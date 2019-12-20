Requirements
============

   This installation requires the pre-installation of the cmake tool,
   a C++11 ready compiler such as g++ version 4.7, and the library: sdsl.

   For Linux, the library sdsl will be installed while makin the tool itself as given below.


Basic Instructions
==================

   To compile the tool, use the following shell command:
```
mkdir build
cd build
cmake ..
make
```
## Usage of the tool: 
 Usage: bin/fraco <args>
 Mandatory args:
  -i, --input-file <str> 	 Input file  name containing fragments (see README.md for the format).
  -o, --output-file <str> 	 Output filename (see README.md for the format).

 Optional args:
  -t, --threads <int> 	 	 Number of threads (DEFAULT: 1). 
  -m, --match <int> 	 	 Score for matching bases (DEFAULT: 5).
  -x, --mismatch <int> 	 	 Score for matching bases (DEFAULT: -4).
  -g, --gap <int> 	 	 Gap penalty (must be negative) (DEFAULT: -8).
  -h, --help 	 	 	 Prints the usage.

 **Example:**  bin/fraco -i ../sample/NC_000913.anc -o ../sample/sample_output.cons

