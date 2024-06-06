# CONSULT-classify

## About CONSULT-classify
CONSULT-classify is a taxonomic classifer, and developed as an extension of [CONSULT-II](https://github.com/bo1929/CONSULT-II).
It works with CONSULT-II's output, which contains a list of pairs for each read. Each pair is distance and taxomic ID, describing a *k*-mer matches of that read.
For example, below is an instace of such an entry for read `@DHUI01000028.1-1174055`.
```
@DHUI01000028.1-1174055
  -- 1671650:2 1671650:2 1671650:4 1671650:2 1671650:2
  rc 1671650:1 1671650:2 1671650:2 1263026:4 1671650:3 1671650:0
```
Here, `rc` stands for reverse-complement form and `--` stands for original form, `1671650:2` is a taxonomic ID and distance pair, respectively.

CONSULT-classify has four argument `--input-matches-dir` (`-i`), `--output-predictions-dir` (`-o`), `--taxonomy-path` and `--thread-count`.
See example experiments to get started.

## Installation
In order to compile the code, `tbb` library must be installed and builtfirst.
Run `misc/get_tbb.sh` to install and build the `tbb` library into the repository's main directory (`lib/tbb`).
You may need to check the paths given for `TBB_LIBRARY_RELEASE` and `TBB_LIBRARY_DEBUG` in the makefile after building `tbb`.
If they are not matching with the build directories in `lib/tbb/build/`, change them accordingly.
Then run `make all` in the main reposity directory to compile the code to genrate the executible `consult_classify`.

## Running experiments
To run `consult_classify` with a benchmark dataset located in `experiments/data`, the `makefile` in the `experiments` directory can be used.
Simply change current directory to `experiments`, and run `make`.
To run sequential version of the `consult_classify`, use `make sequential` and similarly, run `make parallel` to benefit from parallelism.
`make parallel` displays a prompt and takes user input to determine total number of threads, enter the thread count that you want to use and press enter.
Results, i.e., predicitions for the reads in the dataset will be stored in `experiments/predicitions`.
The first column is read-ID, second columns stands for being reverse-complement or not, third column is the predicted NCBI taxonomic ID, and the other two columns are assigned taxonomic ID's vote value & total vote value, respectively.
The total run-times can be found in the log files located in the `experiments/log` directory.
