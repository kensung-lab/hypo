# Suk: Solid (Unique) Kmers Module
=======================================================

Module for finding Solid (Unique) Kmers from genomic reads-data. 

## Dependencies
- Either Mac OS X or Linux are currently supported.
- gcc 4.8+ or clang 3.4+
- cmake 3.8+
- [KMC3](https://github.com/refresh-bio/KMC)
  + Can easily be installed using conda as follows: 
  ```console
  conda install -c bioconda kmc
  ``` 
  + KMC should be in the path (`$PATH`) for Suk to work.

## Building
CmakeLists is provided in the project root folder.
### Building Library
Run the following commands to build a library `suk.a` in `build/lib` :
```console
  git clone --recursive https://github.com/Ritu-Kundu/suk suk
  cd suk
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
  make
```
### Building Executable
Run the following commands to build a binary (executable) `suk` in `build/bin` :
```console
  git clone --recursive https://github.com/Ritu-Kundu/suk suk
  cd suk
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Dsuk_build_executable=ON -Doptimise_for_native=ON ..
  make
```
**Notes:** 
* If `--recursive` was omitted from `git clone`, please run `git submodule init` and `git submodule update` before compiling.
* If target machine is different from the one on which Hypo is being compiled, exclude the flag `-Doptimise_for_native=ON`.

## Usage of the Tool
 ```
 suk <args>

Mandatory args:
  -i, --input <str>   Input file  name containing reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.
  -k, --kmer-len <int>    The length of the kmer (must be positive).
 Optional args:
  -o, --output-prefix <str>   Prefix of the output filename. [DEAFULT: SUK_k<x> where x is the value of k] 
  -t, --threads <int>   Number of threads. [DEFAULT: 1] 
  -d, --dump-txt    Dumps (in addition to the bit-vector) the unique kmers into textfile (<prefix>.txt where prefix is the value of output-prefix). [DEFAULT: False].
  -e, --exclude-hp    Excludes the kmers having homopolymer terminals (see README.md). [DEFAULT: False].
  -c, --coverage <int> 		 Expected coverage of the reads. KMC counter will be set to 4 x expected coverage. [DEFAULT: 50].
  -m, --kmc-memory <int> 	 Max memory usage of KMC (in GB). [DEFAULT: 12].

  -h, --help    Prints the usage.

  ```
### Output Files
- <prefix>.bv file will be created that stores the bit-vector representing the unique kmers. If no `--output-prefix` (or `-o`) is provided, `SUK_k<x>.bv` will be created where `<x>` is the provided value of `k`.
- If text-dump has been activated (using `-d` or `--dump-txt`), <prefix>.txt file will also be created that contains the unique kmers (one in each line). 
 
### Example  
The tool can be run on the sample data (EColi K12) supplied in the `sample` folder, as follows (assuming you are in `build` folder):
 ```console
  ./bin/suk -k 15 -i @../sample/illumina_files.txt
  ```

Here, there are two files for Illumina paired-end reads: SRR1770413_1.fastq and SRR1770413_2.fastq. The names of these two files are written in `illumina_files.txt`. The tool will generate and store the bit-vector representing unique kmers of length 15 in the file `SUK_k15.bv` in the folder from which the tool has been run.


To dump the kmers in text-format, run the following:
```console
  ./bin/suk -d -k 15 -i @../sample/illumina_files.txt
```
This will generate an additional file `SUK_k15.txt` in the folder from which the tool has been run.


## Usage of the Library
You can use Suk in your project by linking `suk` library (See `Building Library (above)`).
An example usage of the library is as follows:
```console
#include "suk/SolidKmers.hpp"

int main(int argc, char **argv) {

  uint32_t threads = 1;
  unsigned int k = 7;
  uint32_t max_memory = 12;
  uint32_t coverage = 50;
  std::vector<std::string> read_filenames;
  read_filenames.push_back("file1"); // first file containing reads

  suk::SolidKmers sk(k);  
  sk.initialise(read_filenames,threads,max_memory, coverage, false);

  // Store the bit-vector into a file 
  sk.store("suk.bv");

  // Load the bit-vector from some previously stored file
  // suk::SolidKmers sk2(k);
  // sk2.load("pvs_stored_file.bv");  

  // Output the solid unique kmers in text format 
  sk.dump_txt("suk.txt");

  /* Check if an encoded kmer is a solid unique kmer */
  uint64_t encoded_kmer = 4u;
  if (sk.is_solid(encoded_kmer)) {std::cout << "Kmer is present\n";}

  /* Check if a kmer (in text format) is a solid unique kmer */
  std::string kmer("AACGTTT"); // kmer of length 7 (=k)
  if (sk.is_solid(kmer)) {std::cout << "Kmer is present\n";}
  

  return 0;
}
```
### Add Suk into your project via CMake
If you are using CMake for compilation, we recommend adding `suk` as a git submodule with the command `git submodule add https://github.com/Ritu-Kundu/suk external/suk`. Afterwards, modify your top level `CMakeLists.txt` file accordingly:
```console
add_subdirectory(external/suk EXCLUDE_FROM_ALL)
target_link_libraries(your_exe suk)
```
Add `#include "suk/SolidKmers.hpp"` to access the suk API in your source file.

## Description
### Kmer Representation
A kmer is a substring of length equal to k. DNA alphabet has four letters/bases:"ACGT". Therefore, there are `4^k` number of distinct kmers possible. DNA alphabet can be represented by 2 bits per base: `00`, `01`, `10`, `11`. Thus each kmer can be represented in `2k` bits. If we encode each base of a kmer into bits, it represents a unique number. Consequently, we have one-to-one correspondence between kmers and numbers of `2k` bits. For example, let k=2. 
```
k=2
# distinct kmers: 4^k = 4^2 = 16   
distinct numbers of 2k bits: 2^(2k) = 2^4 = 16
kmers:          AA    AC    AG    AT    CA    CC    CG    CT    GA    GC    GG    GT    TA    TC    TG    TT  
Encoded kmers:  0000  0001  0010  0011  0100  0101  0110  0111  1000  1001  1010  1011  1100  1101  1110  11111
Corres Number:  0     1     2     3     4     5     6     7     8     9     10    11    12    13    14    15
```
Thus, we can show the collection of kmers appearing in data as a bit-vector of length `2k`. Each index of the vector represents a kmer; if a kmer is present, we set the bit at the corresponding index; otherwise, reset.

### Solid (Unique) kmers
The solid (unique) kmers are those that are expected to appear exactly once in the genome from which reads have sequenced. These can be determined with the help of frequency distribution (histogram) of kmers in the reads. More about the kmer-distribution of a typical genome can be read [here](https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/#).

### Kmer with Homopolymer Terminals
A **homopolymer** is a sequence consisting of the same letter repeated multiple times (more than once), e.g. "AAA". Here, a kmer __starts__ in a homopolymer if the first base is same as the second base. Similarly, a kmer __ends__ in a homopolymer if the last base is same as the second-last base. If exclude-homopolymer filter is on (using `--exclude-hp` or `-e`), we ignore the kmers that start or end in homopolymers.

### Algorithm Steps
- The tool produces the frequency histogram (using [KMC3](https://github.com/refresh-bio/KMC)). 
  + For only the kmers with frequency in the range [2,4*coverage]
- It then finds out the threshold frequencies ([lower,upper]) from the histogram so that the kmers with frequency appearing in the range between the thresholds (inclusive) are expected to be unique.
  + It first finds out the threshold (error-threshold) such that the kmers with frequency below it are most probably the ones containing errors. Data with frequency below error-threshold is discarded from further consideration.
  + Next, it finds out the global-maxima. The frequency at the global-maxima will be taken as the mean-coverage.
  + Next, it computes the lower cut-off by scanning on the left of maxima in a look-forward sliding window fashion: Lower cut-off will be a frequency where the number of kmers is lower than most of those with frequency in an adjacent window (of size 5).
  + Finally, the upper threshold is determined in the same fashion as that of the lower-threshold but on the right-side of the global maxima.
    * There may be some cases where the upper threshold could not be determined in this fashion (because of no uphill curve after maxma); Plan B for such cases is to find the first platue; it works as follows: For every point take average of change (%) from the current value in a window of next 5 values (considering only decreased values). Run moving average (considering a window of size 5) on this range and take the frequency corresponding to the first minimum moving average.
- Lastly, it sets the bits corresponding to all kmers (computed by ntHits) that have frequency <=upper-threshold.
  + Before setting, it checks whether the kmer has HomoPolymer ends (the first and second bases are the same or the last and the second last bases are the same). It does not include kmer which has Homopolymer ends, if exclude-homopolymer filter is on (using `--exclude-hp` or `-e`).

### Notes
The tool
- Assumes the data to be consisting of DNA alphabet: `A`,`C`,`G`,`T`,`N`; (Case insensitive)
- Ignores kmers containing `N`. 
- Kmer-computation is in canonical form (lexicographically smaller amongst the kmer and its reverse-complement). But in the bit-vector, both the solid (unique) kmer and its reverse-complement are set.



## External Libraries

 * [sdsl](https://github.com/simongog/sdsl-lite) library has been used for bit-vector implementation.
 * [slog](https://github.com/Ritu-Kundu/slog) has been used to print time and memory usage.



