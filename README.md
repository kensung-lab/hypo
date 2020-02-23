# HyPo: Super Fast & Accurate Polisher for Long Read Genome Assemblies
=======================================================================

HyPo--a **Hy**brid **Po**lisher-- utilises short as well as long reads within a single run to polish a long reads assembly of small and large genomes. It exploits unique genomic kmers to selectively polish segments of contigs using partial order alignment of selective read-segments. As demonstrated on human genome assemblies, Hypo generates significantly more accurate polished assembly in about one-third time with about half the memory requirements in comparison to contemporary widely used polishers like Racon. 

Please note that "short reads" doesn't necessarily have to be NGS short reads; HiFi genomic reads (e.g. CCS) like those generated from PacBio SequelII could also be used instead. The requirement is that those reads should be highly accurate (>98% accuracy).

Hypo requires the following as input:
+ Short reads/HiFi reads (in FASTA/FASTQ format; can be compressed)
+ Draft contigs (in FASTA/FASTQ format; can be compressed)
+ Alignments between short reads (or HiFi reads) and the draft (in sam/bam format; should contain CIGAR)
  - If Long (noisy) reads are also to be used for polishing, then alignments between long reads and the draft (in sam/bam format; should contain CIGAR)
+ Expected mean coverage of short reads (or HiFi reads) and approximate size of the genome.

In what follows, short reads can be replaced with HiFi reads.

Broadly, we (conceptually) divide a draft (uncorrected) contig into two types of regions (segments): *Strong* and *Weak*. Strong regions are those which have strong evidence (*support*) of their correctness and thus do not need polishing. Weak regions, on the other hand, will be polished using POA. Each weak region will be polished using either short reads or long reads; short reads taking precedence over long reads. To identify strong regions, we make use of *solid* kmers (expected unique genomic kmers). Strong regions also play a role in selecting the read-segments to polish their neighbouring weak regions. Furthermore, our approach takes into account that the long reads and thus the assemblies generated from them are prone to homopolymer errors as mentioned in the beginning.

## Installation
Hypo is only available for Unix-like platforms (Linux and MAC OS). We recommend using *Option_1* for a convenient installation and in the case where target machine is different from the one on which hypo is installed. On the other hand, *Option_2* is more suitable if the binary of hypo is to be run on the same machine on which it is compiled because then a machine-specific optimised (and thus slightly faster) binary can be produced using the flag `-Doptimise_for_native=ON`.


### Option_1: Conda Package Installation
The convenient way of installation is using the conda package as follows:
```console
conda install -c bioconda hypo
```
Htslib 1.10 may cause conflicts with some of the already installed packages. In such a case, hypo can be installed in a new environment as follows:
```console
conda create --name hypo_env
conda activate hypo_env  
conda install -c bioconda hypo
```

### Option_2: Installation from the source
CmakeLists is provided in the project root folder. 

#### Pre-requisites
For installing from the source, the following requirements are assumed to be installed already (with path to their binaries available in $PATH).
- Zlib
- OpenMP
- GCC (>=7.3)
  * Following are the commands to update GCC (say to GCC 8) on an Ubuntu machine (from say GCC 5):
  ```console
    sudo apt-get update; sudo apt-get install build-essential software-properties-common -y;
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y; sudo apt update; 
    sudo apt install gcc-snapshot -y; sudo apt update
    sudo apt install gcc-8 g++-8 -y; 
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 60 --slave /usr/bin/g++ g++ /usr/bin/g++-8
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5;
  ```
- [HTSLIB](https://github.com/samtools/htslib) (version >=1.10)
  + If htslib version 1.10 or higher is not installed, we recommend using `install_deps.sh` in the project folder to install it locally.

- [KMC3](https://github.com/refresh-bio/KMC)
  + Required for suk
  + Can also be installed using conda as follows: 
  ```console
  conda install -c bioconda kmc
  ``` 

#### Building Executable
Run the following commands to build a binary (executable) `hypo` in `build/bin` :
If htslib version 1.10 or higher is installed:
```console
  git clone --recursive https://github.com/kensung-lab/hypo hypo
  cd hypo
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
  make -j 8
```
If htslib is not installed or the version is smaller than 1.10:
```console
  git clone --recursive https://github.com/kensung-lab/hypo hypo
  cd hypo
  chmod +x install_deps.sh
  ./install_deps.sh
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
  make -j 8
```
**Notes:** 
* If `--recursive` was omitted from `git clone`, please run `git submodule init` and `git submodule update` before compiling.
* If target machine is different from the one on which Hypo is being compiled, exclude the flag `-Doptimise_for_native=ON`.


## Usage of the tool: 
```console
 Usage: hypo <args>

 ** Mandatory args:
	-r, --reads-short <str>
	Input file name containing reads (in fasta/fastq format; can be compressed). A list of files containing file names in each line can be passed with @ prefix.

	-d, --draft <str>
	Input file name containing the draft contigs (in fasta/fastq format; can be compressed). 

	-b, --bam-sr <str>
	Input file name containing the alignments of short reads against the draft (in bam/sam format; must have CIGAR information). 

	-c, --coverage-short <int>
	Approximate mean coverage of the short reads. 

	-s, --size-ref <str>
	Approximate size of the genome (a number; could be followed by units k/m/g; e.g. 10m, 2.3g). 

** Optional args:
	-B, --bam-lr <str>
	Input file name containing the alignments of long reads against the draft (in bam/sam format; must have CIGAR information). 
	[Only Short reads polishing will be performed if this argument is not given]

	-o, --output <str>
	Output file name. 
	[Default] hypo_<draft_file_name>.fasta in the working directory.

 	-t, --threads <int>
	Number of threads. 
	[Default] 1.

 	-p, --processing-size <int>
	Number of contigs to be processed in one batch. Lower value means less memory usage but slower speed. 
	[Default] All the contigs in the draft.

	-k, --kind-sr <str>
	Kind of the short reads. 
	[Valid Values] 
		sr	(Corresponding to NGS reads like Illumina reads) 
		ccs	(Corresponding to HiFi reads like PacBio CCS reads) 
	[Default] sr.


 	-m, --match-sr <int>
	Score for matching bases for short reads. 
	[Default] 5.

 	-x, --mismatch-sr <int>
	Score for mismatching bases for short reads. 
	[Default] -4.

 	-g, --gap-sr <int>
	Gap penalty for short reads (must be negative). 
	[Default] -8.

 	-M, --match-lr <int>
	Score for matching bases for long reads. 
	[Default] 3.

 	-X, --mismatch-lr <int>
	Score for mismatching bases for long reads. 
	[Default] -5.

 	-G, --gap-lr <int>
	Gap penalty for long reads (must be negative). 
	[Default] -4.

 	-n, --ned-th <int>
	Threshold for Normalised Edit Distance of long arms allowed in a window (in %). Higher number means more arms allowed which may slow down the execution.
	[Default] 20.

 	-q, --qual-map-th <int>
	Threshold for mapping quality of reads. The reads with mapping quality below this threshold will not be taken into consideration. 
	[Default] 2.

 	-i, --intermed
	Store or use (if already exist) the intermediate files. 
	[Currently, only Solid kmers are stored as an intermediate file.].

 	-h, --help
	Print the usage. 


```

### Output File
If no `--output` (or `-o`) is provided, `hypo_X.fasta` will be created in the working directory where <X> is the name of the draft file. 

### Intermediate Files
If `--intermed` (or `-i`) is used, hypo will store the intermediate files corresponding to the solid kmers in a folder named `aux` in the first run. Another file indicating the progress of the run will also be created in that folder. In the next run (with `-i`), those intermediate files will be used instead of running `suk` module. Remove `-i` or delete `aux` folder to start Hypo from the beginning. **It is recommended to use this flag while going for multiple rounds of polishing using Hypo so as to make the subsequent rounds faster.**

### Resource Usage
The option `--processing-size` (or `-p`) controls the number of contigs processed in one batch. By default, all contigs will be processed in a single batch. More the number of contigs processed in a batch higher will be the memory used. Lower number, on the other hand, may not utilise the number of threads efficiently. As a reference, for the whole human genome we used `-p 96` on a 48 core machine which used about 380G RAM and finished its run in about 3 hours (only Illumina polishing). We recommend using `-p` to be a multiple of `-t`. If the genome size is not too large for the available RAM, we recommend processing all the contigs in a single batch (i.e. avoid specifying `-p`).

### Normalised Edit Distance
The option `--ned-th` (or `-n`) is relevant only when long reads polishing is also being done. It sets the threshold for the Normalised Edit Distance which measures similarity between a long read and the corresponding draft segment to which it is aligned. Formally, Normalised edit distance = (1bp edit distance between aligned read_segment and draft /aligned_length of draft)x100
Higher the threshold, more read-segments will be taken into consideration. We recommend using the default value for the whole human genome polishing. 
 
### Example 1 (using Illumina paired-end reads)

#### Mapping the short reads to contigs:
Assuming `$R1`, `$R2`, `$Draft` contain the names of the files containing short reads (paired-end) and draft contigs, respectively. Let `$NUMTH` represents the number of threads  to be used.
```console
minimap2 --secondary=no --MD -ax sr -t $NUMTH $DRAFT $R1 $R2 | samtools view -Sb - > mapped-sr.bam
samtools sort -@$NUMTH -o mapped-sr.sorted.bam mapped-sr.bam
samtools index mapped-sr.sorted.bam
rm mapped-sr.bam
```

#### Mapping the long reads to contigs:
Assuming `$LONGR` and `$Draft` contain the names of the files containing the long reads (PacBio or ONT) and draft contigs, respectively. Let `$NUMTH` represents the number of threads to be used. `$RTYPE` will be `map-pb` for PacBio and `map-ont` for ONT reads.
```console
minimap2 --secondary=no --MD -ax $RTYPE -t $NUMTH $DRAFT $LONGR | samtools view -Sb - > mapped-lg.bam
samtools sort -@$NUMTH -o mapped-lg.sorted.bam mapped-lg.bam
samtools index mapped-lg.sorted.bam
rm mapped-lg.bam
```

#### Running Hypo:
Create a text file containing the names of the short reads files.
```console
echo -e "$R1\n$R2" > il_names.txt
```

Let genome size be around 3Gbp and average coverage of Illumina reads is around 55. Say, we would like to use 48 threads and process 96 contigs in a single batch.

A sample run of Hypo (for only short reads polishing) can be:
```console
./bin/hypo -d $DRAFT -r @il_names.txt -s 3g -c 55 -b mapped-sr.sorted.bam -p 96 -t 48 -o whole_genome.h.fa
```

A sample run of Hypo (for short reads as well as long reads polishing) can be:
```console
./bin/hypo -d $DRAFT -r @il_names.txt -s 3g -c 55 -b mapped-sr.sorted.bam -B mapped-lg.sorted.bam -p 96 -t 48 -o whole_genome.h2.fa
```
### Example 2 (using CCS reads)

#### Mapping the CCS reads to contigs:
Assuming `$READS` and `$Draft` contain the names of the files containing CCS reads and draft contigs, respectively. Let `$NUMTH` represents the number of threads  to be used.
```console
minimap2 --secondary=no --MD -ax asm20 -t $NUMTH $DRAFT $READS | samtools view -Sb - > mapped-ccs.bam
samtools sort -@$NUMTH -o mapped-ccs.sorted.bam mapped-ccs.bam
samtools index mapped-ccs.sorted.bam
rm mapped-ccs.bam
```

#### Running Hypo:
Let genome size be around 3Gbp and average coverage of CCS reads is around 30. Say, we would like to use 48 threads and process all the contigs in a single batch.
A sample run of Hypo (for only CCS polishing) can be:
```console
./bin/hypo -d $DRAFT -r $READS -s 3g -c 30 -b mapped-ccs.sorted.bam -t 48 -o whole_genome.h.fa
```

## Method and Results
For the whole human genome (HG002) assembly, Hypo took about 3 hours and about 380G RAM to polish (on a 48 cores machine with 48 threads) using only Illumina reads. For polishing using Illumina as well as PacBio reads, time taken was about 4 hours and 15 minutes using about 410G RAM. The method and partial results can be found at [BioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.19.882506v1).
```
HyPo: Super Fast & Accurate Polisher for Long Read Genome Assemblies
Ritu Kundu, Joshua Casey, Wing-Kin Sung
bioRxiv 2019.12.19.882506; doi: https://doi.org/10.1101/2019.12.19.882506 
```

## Limitations
 * N's in the contig/scaffold to be polished are considered as any other base (i.e. an *N* will be replaced by the most frequent base aligned to it in the corresponding POA graph). 
 * The maximum length a contig/scaffold can have is limited to 4,294,967,295 (the maximum value for a 32-bit unsigned integer). 


## External Libraries

 * [sdsl-lite](https://github.com/simongog/sdsl-lite) has been used for rank-select and bit-vectors data-structures.
 * [slog](https://github.com/Ritu-Kundu/slog) has been used to print time and memory usage.
 * [suk](https://github.com/Ritu-Kundu/suk) has been used as the module to compute the solid (unique) kmers.
 * An adapted version of [spoa](https://github.com/rvaser/spoa.git) library has been used for POA.

## Contact
Other than raising issues on Github, you can contact Ritu Kundu (dcsritu@nus.edu.sg) or Joshua Casey (joshuac@comp.nus.edu.sg) for getting help in installation/usage or any other related query.





