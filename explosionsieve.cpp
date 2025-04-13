#include "explosionsieve.hpp"

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "external/chess.hpp"
#include "external/gzip/gzstream.h"
#include "external/parallel_hashmap/phmap.h"
#include "external/threadpool.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

using namespace chess;

// we want to store the ratio and the move comment, plus other info
using Statistics = std::pair<float, std::string>;

using fen_map_t = phmap::parallel_flat_hash_map<
    std::string, Statistics, std::hash<std::string>, std::equal_to<std::string>,
    std::allocator<std::pair<const std::string, Statistics>>, 8, std::mutex>;

// map to hold move counters that cutechess-cli changed from original FENs
using map_fens = std::unordered_map<std::string, std::pair<int, int>>;

fen_map_t fen_map;

// map to collect metadata for tests
using map_meta = std::unordered_map<std::string, TestMetaData>;

std::atomic<std::size_t> total_files = 0;
std::atomic<std::size_t> total_games = 0;
std::atomic<float> min_ratio = std::numeric_limits<float>::infinity();
;

namespace analysis {

class Analyze : public pgn::Visitor {
public:
  Analyze(std::string_view file, const map_meta &meta_map,
          const std::string &regex_engine, const map_fens &fixfen_map,
          float ratio_bound, float bf, std::mutex &progress_output)
      : file(file), meta_map(meta_map), regex_engine(regex_engine),
        fixfen_map(fixfen_map), ratio_bound(ratio_bound), bf(bf),
        progress_output(progress_output) {}

  virtual ~Analyze() {}

  void startPgn() override {}

  void header(std::string_view key, std::string_view value) override {

    if (key == "Termination") {
      if (value == "time forfeit" || value == "abandoned" ||
          value == "stalled connection" || value == "illegal move" ||
          value == "unterminated") {
        valid_game = false;
      }
    }

    if (key == "FEN") {
      std::regex p("^(.+) 0 1$");
      std::smatch match;
      std::string value_str(value);

      // revert changes by cutechess-cli to move counters
      if (!fixfen_map.empty() && std::regex_search(value_str, match, p) &&
          match.size() > 1) {
        std::string fen = match[1];
        auto it = fixfen_map.find(fen);

        if (it == fixfen_map.end()) {
          std::cerr << "While parsing " << file << " could not find FEN " << fen
                    << " in fixFENsource." << std::endl;
          board.setFen(value);
        } else {
          const auto &fix = it->second;
          std::string fixed_value = fen + " " + std::to_string(fix.first) +
                                    " " + std::to_string(fix.second);
          board.setFen(fixed_value);
        }
      } else
        board.setFen(value);
    }

    if (key == "Variant" && value == "fischerandom") {
      board.set960(true);
    }

    if (key == "White") {
      white = value;
    }

    if (key == "Black") {
      black = value;
    }

    if (key == "WhiteElo") {
      whiteElo = std::atoi(value.data());
    }
    if (key == "BlackElo") {
      blackElo = std::atoi(value.data());
    }

    if (key == "TimeControl") {
      fs::path path(file);
      std::string filename = path.filename().string();
      std::string test_id = filename.substr(0, filename.find_first_of("-."));
      std::string test_filename = (path.parent_path() / test_id).string();
      if (meta_map.find(test_filename) != meta_map.end() &&
          meta_map.at(test_filename).tc.has_value()) {

        std::string meta_str(meta_map.at(test_filename).tc.value());
        const size_t metaslash_pos = meta_str.find("/");
        const size_t metastart_pos =
            metaslash_pos == std::string::npos ? 0 : metaslash_pos + 1;
        const size_t metaplus_pos = meta_str.find("+");
        const auto metamatch_basetc =
            meta_str.substr(metastart_pos, metaplus_pos);
        float metabasetc = fast_stof(metamatch_basetc.data());

        std::string value_str(value);
        const size_t slash_pos = value_str.find("/");
        const size_t start_pos =
            slash_pos == std::string::npos ? 0 : slash_pos + 1;
        const size_t plus_pos = value_str.find("+");
        const auto match_basetc = value_str.substr(start_pos, plus_pos);
        float basetc = fast_stof(match_basetc.data());

        time_control = "'" + value_str + "' '" + meta_str + "'";
        tc_factor = metabasetc / basetc;
      } else {
        std::cerr << "While parsing " << file
                  << " could not find meta data with TC." << std::endl;
      }
    }
  }

  void startMoves() override {
    if (tc_factor == 0.0 || !valid_game) {
      this->skipPgn(true);
      return;
    }
    do_filter = !regex_engine.empty();

    if (do_filter) {
      if (white.empty() || black.empty()) {
        this->skipPgn(true);
        return;
      }

      std::regex regex(regex_engine);

      if (std::regex_match(white, regex)) {
        filter_side = Color::WHITE;
      }

      if (std::regex_match(black, regex)) {
        if (filter_side == Color::NONE) {
          filter_side = Color::BLACK;
        } else {
          do_filter = false;
        }
      }
    }
    total_games++;
  }

  void move(std::string_view move, std::string_view comment) override {

    // fishtest uses Nf3 {+0.57/17 2.313s}
    const size_t slash_pos = comment.find("/");
    const size_t space_pos = comment.find(" ");

    if (!do_filter || filter_side == board.sideToMove()) {
      if (comment != "book" && slash_pos != std::string::npos &&
          space_pos != std::string::npos && slash_pos < space_pos) {

        const auto match_depth = comment.substr(slash_pos + 1, space_pos);
        const auto match_time =
            comment.substr(space_pos + 1, comment.size() - 1);
        int depth = std::max(1, std::stoi(match_depth.data()));
        float time = std::max(float(1e-16), fast_stof(match_time.data()));

        float ratio =
            (bf == 0.0 ? depth : std::pow(depth, bf)) / (time * tc_factor);
        if (ratio < ratio_bound) {
          auto key = board.getFen();
          auto engine = board.sideToMove() == Color::WHITE ? white : black;
          auto str = "'" + std::string(comment) + "' " + time_control + " '" +
                     engine + "'";
          auto val = std::pair<float, std::string>(ratio, str);
          fen_map.lazy_emplace_l(
              std::move(key),
              [&](fen_map_t::value_type &p) {
                if (ratio < p.second.first) {
                  p.second = val;
                  if (ratio < min_ratio)
                    min_ratio = ratio;
                }
              },
              [&](const fen_map_t::constructor &ctor) {
                ctor(std::move(key), val);
                if (ratio < min_ratio)
                  min_ratio = ratio;
              });
        }
      }
    }

    try {
      Move m = uci::parseSan(board, move, moves);

      // chess-lib may call move() with empty strings for move
      if (m == Move::NO_MOVE) {
        this->skipPgn(true);
        return;
      }

      board.makeMove<true>(m);
    } catch (const uci::AmbiguousMoveError &e) {
      std::cerr << "While parsing " << file << " encountered: " << e.what()
                << '\n';
      this->skipPgn(true);
    }
  }

  void endPgn() override {
    board.set960(false);
    board.setFen(constants::STARTPOS);

    filter_side = Color::NONE;

    white.clear();
    black.clear();
    time_control.clear();

    tc_factor = whiteElo = blackElo = 0;
    valid_game = true;
  }

private:
  std::string_view file;
  const map_meta &meta_map;
  const std::string &regex_engine;
  const map_fens &fixfen_map;
  float ratio_bound, bf;
  std::mutex &progress_output;

  Board board;
  Movelist moves;

  bool skip = false;

  bool do_filter = false;
  Color filter_side = Color::NONE;

  std::string white, black, time_control;

  int whiteElo = 0, blackElo = 0;
  float tc_factor = 0.0;
  bool valid_game = true;
};

void ana_files(const std::vector<std::string> &files, const map_meta &meta_map,
               const std::string &regex_engine, const map_fens &fixfen_map,
               float ratio_bound, float bf, std::mutex &progress_output) {

  for (const auto &file : files) {
    std::string move_counter;
    const auto pgn_iterator = [&](std::istream &iss) {
      auto vis =
          std::make_unique<Analyze>(file, meta_map, regex_engine, fixfen_map,
                                    ratio_bound, bf, progress_output);

      pgn::StreamParser parser(iss);

      try {
        parser.readGames(*vis);
      } catch (const std::exception &e) {
        std::cout << "Error when parsing: " << file << std::endl;
        std::cerr << e.what() << '\n';
      }
    };

    if (file.size() >= 3 && file.substr(file.size() - 3) == ".gz") {
      igzstream input(file.c_str());
      pgn_iterator(input);
    } else {
      std::ifstream pgn_stream(file);
      pgn_iterator(pgn_stream);
      pgn_stream.close();
    }

    ++total_files;

    // Limit the scope of the lock
    {
      const std::lock_guard<std::mutex> lock(progress_output);
      std::cout << "\rProcessed " << total_files << " files, found "
                << fen_map.size()
                << " positions with lowest ratio = " << min_ratio << std::flush;
    }
  }
}

} // namespace analysis

[[nodiscard]] map_fens get_fixfen(std::string file) {
  map_fens fixfen_map;
  if (file.empty()) {
    return fixfen_map;
  }

  const auto fen_iterator = [&](std::istream &iss) {
    std::string line;
    while (std::getline(iss, line)) {
      std::istringstream iss(line);
      std::string f1, f2, f3, ep;
      int halfmove, fullmove = 0;

      iss >> f1 >> f2 >> f3 >> ep >> halfmove >> fullmove;

      if (!fullmove)
        continue;

      auto key = f1 + ' ' + f2 + ' ' + f3 + ' ' + ep;
      auto fixfen_data = std::pair<int, int>(halfmove, fullmove);

      if (fixfen_map.find(key) != fixfen_map.end()) {
        // for duplicate FENs, prefer the one with lower full move counter
        if (fullmove < fixfen_map[key].second) {
          fixfen_map[key] = fixfen_data;
        }
      } else {
        fixfen_map[key] = fixfen_data;
      }
    }
  };

  if (file.size() >= 3 && file.substr(file.size() - 3) == ".gz") {
    igzstream input(file.c_str());
    fen_iterator(input);
  } else {
    std::ifstream input(file);
    fen_iterator(input);
  }

  return fixfen_map;
}

[[nodiscard]] map_meta get_metadata(const std::vector<std::string> &file_list,
                                    bool allow_duplicates) {
  map_meta meta_map;
  // map to check for duplicate tests
  std::unordered_map<std::string, std::string> test_map;
  std::set<std::string> test_warned;

  for (const auto &pathname : file_list) {
    fs::path path(pathname);
    std::string filename = path.filename().string();
    std::string test_id = filename.substr(0, filename.find_first_of("-."));
    std::string test_filename = (path.parent_path() / test_id).string();

    if (test_map.find(test_id) == test_map.end()) {
      test_map[test_id] = test_filename;
    } else if (test_map[test_id] != test_filename) {
      if (test_warned.find(test_filename) == test_warned.end()) {
        std::cout << (allow_duplicates ? "Warning" : "Error")
                  << ": Detected a duplicate of test " << test_id
                  << " in directory " << path.parent_path().string()
                  << std::endl;
        test_warned.insert(test_filename);

        if (!allow_duplicates) {
          std::cout << "Use --allowDuplicates to continue nonetheless."
                    << std::endl;
          std::exit(1);
        }
      }
    }

    // load the JSON data from disk, only once for each test
    if (meta_map.find(test_filename) == meta_map.end()) {
      std::ifstream json_file(test_filename + ".json");

      if (!json_file.is_open())
        continue;

      json metadata = json::parse(json_file);

      meta_map[test_filename] = metadata.get<TestMetaData>();
    }
  }

  return meta_map;
}

template <typename STRATEGY>
void filter_files(std::vector<std::string> &file_list, const map_meta &meta_map,
                  const STRATEGY &strategy) {
  const auto applier = [&](const std::string &pathname) {
    fs::path path(pathname);
    std::string filename = path.filename().string();
    std::string test_id = filename.substr(0, filename.find_first_of("-."));
    std::string test_filename = (path.parent_path() / test_id).string();
    return strategy.apply(test_filename, meta_map);
  };
  const auto it = std::remove_if(file_list.begin(), file_list.end(), applier);
  file_list.erase(it, file_list.end());
}

class BookFilterStrategy {
  std::regex regex_book;
  bool invert;

public:
  BookFilterStrategy(const std::regex &rb, bool inv)
      : regex_book(rb), invert(inv) {}

  bool apply(const std::string &filename, const map_meta &meta_map) const {
    // check if metadata and "book" entry exist
    if (meta_map.find(filename) != meta_map.end() &&
        meta_map.at(filename).book.has_value()) {
      bool match =
          std::regex_match(meta_map.at(filename).book.value(), regex_book);
      return invert ? match : !match;
    }

    // missing metadata or "book" entry can never match
    return true;
  }
};

class RevFilterStrategy {
  std::regex regex_rev;

public:
  RevFilterStrategy(const std::regex &rb) : regex_rev(rb) {}

  bool apply(const std::string &filename, const map_meta &meta_map) const {
    if (meta_map.find(filename) == meta_map.end()) {
      return true;
    }

    if (meta_map.at(filename).resolved_base.has_value() &&
        std::regex_match(meta_map.at(filename).resolved_base.value(),
                         regex_rev)) {
      return false;
    }

    if (meta_map.at(filename).resolved_new.has_value() &&
        std::regex_match(meta_map.at(filename).resolved_new.value(),
                         regex_rev)) {
      return false;
    }

    return true;
  }
};

class TcFilterStrategy {
  std::regex regex_tc;

public:
  TcFilterStrategy(const std::regex &rb) : regex_tc(rb) {}

  bool apply(const std::string &filename, const map_meta &meta_map) const {
    if (meta_map.find(filename) == meta_map.end()) {
      return true;
    }

    if (meta_map.at(filename).new_tc.has_value() &&
        meta_map.at(filename).tc.has_value()) {
      if (meta_map.at(filename).new_tc.value() !=
          meta_map.at(filename).tc.value()) {
        return true;
      }

      if (std::regex_match(meta_map.at(filename).tc.value(), regex_tc)) {
        return false;
      }
    }

    return true;
  }
};

class TcEqualFilter {
public:
  TcEqualFilter() {}

  bool apply(const std::string &filename, const map_meta &meta_map) const {
    if (meta_map.find(filename) == meta_map.end()) {
      return true;
    }

    if (meta_map.at(filename).new_tc.has_value() &&
        meta_map.at(filename).tc.has_value()) {
      if (meta_map.at(filename).new_tc.value() !=
          meta_map.at(filename).tc.value()) {
        return true;
      }

      return false;
    }

    return true;
  }
};

class ThreadsFilterStrategy {
  int threads;

public:
  ThreadsFilterStrategy(int t) : threads(t) {}

  bool apply(const std::string &filename, const map_meta &meta_map) const {
    if (meta_map.find(filename) == meta_map.end()) {
      return true;
    }

    if (meta_map.at(filename).threads.has_value() &&
        meta_map.at(filename).threads.value() == threads) {
      return false;
    }

    return true;
  }
};

void process(const std::vector<std::string> &files_pgn,
             const map_meta &meta_map, const std::string &regex_engine,
             const map_fens &fixfen_map, float ratio_bound, float bf,
             int concurrency) {
  // Create more chunks than threads to prevent threads from idling.
  int target_chunks = 4 * concurrency;

  auto files_chunked = split_chunks(files_pgn, target_chunks);

  std::cout << "Found " << files_pgn.size() << " .pgn(.gz) files, creating "
            << files_chunked.size() << " chunks for processing." << std::endl;

  // Mutex for progress output
  std::mutex progress_output;

  // Create a thread pool
  ThreadPool pool(concurrency);

  for (const auto &files : files_chunked) {

    pool.enqueue([&files, &meta_map, &regex_engine, &fixfen_map, &ratio_bound,
                  &bf, &progress_output, &files_chunked]() {
      analysis::ana_files(files, meta_map, regex_engine, fixfen_map,
                          ratio_bound, bf, progress_output);
    });
  }

  // Wait for all threads to finish
  pool.wait();
}

void print_usage(char const *program_name) {
  std::stringstream ss;

  // clang-format off
    ss << "Usage: " << program_name << " [options]" << "\n";
    ss << "Options:" << "\n";
    ss << "  --file <path>          Path to .pgn(.gz) file" << "\n";
    ss << "  --dir <path>           Path to directory containing .pgn(.gz) files (default: pgns)" << "\n";
    ss << "  -r                     Search for .pgn(.gz) files recursively in subdirectories" << "\n";
    ss << "  --allowDuplicates      Allow duplicate directories for test pgns" << "\n";
    ss << "  --concurrency <N>      Number of concurrent threads to use (default: maximum)" << "\n";
    ss << "  --matchRev <regex>     Filter data based on revision SHA in metadata" << "\n";
    ss << "  --matchEngine <regex>  Filter data based on engine name in pgns, defaults to matchRev if given" << "\n";
    ss << "  --matchTC <regex>      Filter data based on time control in metadata" << "\n";
    ss << "  --matchThreads <N>     Filter data based on used threads in metadata" << "\n";
    ss << "  --matchBook <regex>    Filter data based on book name" << "\n";
    ss << "  --matchBookInvert      Invert the filter" << "\n";
    ss << "  -o <path>              Path to output epd file (default: explosionsieve.epd)" << "\n";
    ss << "  --fixFENsource         Patch move counters lost by cutechess-cli based on FENs in this file" << "\n";
    ss << "  --ratioBound <rb>      Upper bound for the depth/second ratio (default: 1)" << "\n";
    ss << "  --branchingFactor <bf> If given, replace depth with bf^depth in ratio calculation" << "\n";
    ss << "  --help                 Print this help message" << "\n";
  // clang-format on

  std::cout << ss.str();
}

int main(int argc, char const *argv[]) {
  CommandLine cmd(argc, argv);

  std::vector<std::string> files_pgn;
  std::string default_path = "./pgns";
  std::string regex_engine, regex_book, filename = "explosionsieve.epd";
  map_fens fixfen_map;
  float ratio_bound = 1.0, branching_factor = 0.0;
  std::string branching_factor_str;
  int concurrency = std::max(1, int(std::thread::hardware_concurrency()));

  if (cmd.has_argument("--help", true)) {
    print_usage(argv[0]);
    return 0;
  }

  if (cmd.has_argument("--ratioBound")) {
    ratio_bound = fast_stof(cmd.get_argument("--ratioBound"));
  }
  if (cmd.has_argument("--branchingFactor")) {
    branching_factor_str = cmd.get_argument("--branchingFactor");
    branching_factor = fast_stof(branching_factor_str);
    if (branching_factor <= 0.0) {
      std::cout << "Error: Branching factor must be positive." << std::endl;
      std::exit(1);
    }
  }

  std::string depth_str =
      branching_factor == 0.0 ? "" : ("^" + branching_factor_str);
  std::cout << "Looking for search explosions with depth" << depth_str
            << "/second below " << ratio_bound << " ..." << std::endl;

  if (cmd.has_argument("--concurrency")) {
    concurrency = std::stoi(cmd.get_argument("--concurrency"));
  }

  if (cmd.has_argument("--file")) {
    files_pgn = {cmd.get_argument("--file")};
    if (!fs::exists(files_pgn[0])) {
      std::cout << "Error: File not found: " << files_pgn[0] << std::endl;
      std::exit(1);
    }
  } else {
    auto path = cmd.get_argument("--dir", default_path);

    bool recursive = cmd.has_argument("-r", true);
    std::cout << "Looking " << (recursive ? "(recursively) " : "")
              << "for pgn files in " << path << std::endl;

    files_pgn = get_files(path, recursive);

    // sort to easily check for "duplicate" files, i.e. "foo.pgn.gz" and
    // "foo.pgn"
    std::sort(files_pgn.begin(), files_pgn.end());

    for (size_t i = 1; i < files_pgn.size(); ++i) {
      if (files_pgn[i].find(files_pgn[i - 1]) == 0) {
        std::cout << "Error: \"Duplicate\" files: " << files_pgn[i - 1]
                  << " and " << files_pgn[i] << std::endl;
        std::exit(1);
      }
    }
  }

  std::cout << "Found " << files_pgn.size() << " .pgn(.gz) files in total."
            << std::endl;

  auto meta_map =
      get_metadata(files_pgn, cmd.has_argument("--allowDuplicates", true));

  if (cmd.has_argument("--matchBook")) {
    auto regex_book = cmd.get_argument("--matchBook");

    if (!regex_book.empty()) {
      bool invert = cmd.has_argument("--matchBookInvert", true);
      std::cout << "Filtering pgn files " << (invert ? "not " : "")
                << "matching the book name " << regex_book << std::endl;
      filter_files(files_pgn, meta_map,
                   BookFilterStrategy(std::regex(regex_book), invert));
    }
  }

  if (cmd.has_argument("--matchRev")) {
    auto regex_rev = cmd.get_argument("--matchRev");

    if (!regex_rev.empty()) {
      std::cout << "Filtering pgn files matching revision SHA " << regex_rev
                << std::endl;
      filter_files(files_pgn, meta_map,
                   RevFilterStrategy(std::regex(regex_rev)));
    }

    regex_engine = regex_rev;
  }

  if (cmd.has_argument("--fixFENsource"))
    fixfen_map = get_fixfen(cmd.get_argument("--fixFENsource"));

  if (cmd.has_argument("--matchEngine")) {
    regex_engine = cmd.get_argument("--matchEngine");
  }

  bool applied_tc_filter = false;
  if (cmd.has_argument("--matchTC")) {
    auto regex_tc = cmd.get_argument("--matchTC");

    if (!regex_tc.empty()) {
      std::cout << "Filtering pgn files matching TC " << regex_tc << std::endl;
      filter_files(files_pgn, meta_map, TcFilterStrategy(std::regex(regex_tc)));
      applied_tc_filter = true;
    }
  }

  if (!applied_tc_filter) { // the matchTC filter already checks for this
    std::cout << "Filtering pgn with equal TC." << std::endl;
    filter_files(files_pgn, meta_map, TcEqualFilter());
  }

  if (cmd.has_argument("--matchThreads")) {
    int threads = std::stoi(cmd.get_argument("--matchThreads"));

    std::cout << "Filtering pgn files using threads = " << threads << std::endl;
    filter_files(files_pgn, meta_map, ThreadsFilterStrategy(threads));
  }

  if (cmd.has_argument("-o")) {
    filename = cmd.get_argument("-o");
  }

  std::ofstream out_file(filename);

  const auto t0 = std::chrono::high_resolution_clock::now();

  process(files_pgn, meta_map, regex_engine, fixfen_map, ratio_bound,
          branching_factor, concurrency);

  for (const auto &pair : fen_map) {
    std::string fen = pair.first;
    auto board = Board(fen);
    out_file << fen << " c0 \"d" << depth_str << "/s: " << pair.second.first
             << " " << pair.second.second << "\";\n";
  }
  out_file.close();

  const auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "\nSaved " << fen_map.size() << " unique explosions from "
            << total_games << " games to " << filename << "."
            << "\nTotal time for processing: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                       .count() /
                   1000.0
            << " s" << std::endl;

  return 0;
}
