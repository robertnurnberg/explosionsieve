# Find search explosions in (fishtest) pgn data

Given a set of .pgn(.gz) files with engine evaluations, find the positions
that took the longest to increase search depth.

Example usage:

```
> ./explosionsieve -r --matchTC "60\+0.6" --branchingFactor 1.5 --ratioBound 5
Looking for search explosions with depth^1.5/second below 5 ...
Looking (recursively) for pgn files in ./pgns
Found 865 .pgn(.gz) files in total.
Filtering pgn files matching TC 60\+0.6
Found 799 .pgn(.gz) files, creating 89 chunks for processing.
Processed 799 files, found 1025 positions with lowest ratio = 1.07458
Saved 1025 unique explosions from 80405775 games to explosionsieve.epd.
Total time for processing: 334.188 s
```

```
Usage: ./explosionsieve [options]
Options:
  --file <path>          Path to .pgn(.gz) file
  --dir <path>           Path to directory containing .pgn(.gz) files (default: pgns)
  -r                     Search for .pgn(.gz) files recursively in subdirectories
  --allowDuplicates      Allow duplicate directories for test pgns
  --concurrency <N>      Number of concurrent threads to use (default: maximum)
  --matchRev <regex>     Filter data based on revision SHA in metadata
  --matchEngine <regex>  Filter data based on engine name in pgns, defaults to matchRev if given
  --matchTC <regex>      Filter data based on time control in metadata
  --matchThreads <N>     Filter data based on used threads in metadata
  --matchBook <regex>    Filter data based on book name
  --matchBookInvert      Invert the filter
  -o <path>              Path to output epd file (default: explosionsieve.epd)
  --fixFENsource         Patch move counters lost by cutechess-cli based on FENs in this file
  --ratioBound <rb>      Upper bound for the depth/second ratio (default: 1)
  --branchingFactor <bf> If given, replace depth with bf^depth in ratio calculation
  --help                 Print this help message
```

The code is based on [blundersieve](https://github.com/robertnurnberg/blundersieve).
