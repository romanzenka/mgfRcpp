# mgfRcpp

World's worst Mascot Generic Format parser.

Done to allow loading of very large .mgf files where other packages seemed to fail on out of memory errors / were too slow.

This code has minimal overhead as all the parsing is done using Rcpp, only two resulting data frames are created.
