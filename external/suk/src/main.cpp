/*
 * 
 * Copyright (c) 2019, Ritu Kundu and Joshua Casey
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** Module containing main() method.
 */

#include <getopt.h>
#include <cctype>
#include <math.h> 
#include <algorithm>

#include "suk/SolidKmers.hpp"


namespace suk{

/*
   * Usage of the tool
   */
void usage(void)
{
  std::cout << " Usage: suk <args>\n";
  std::cout << " Mandatory args:\n";
  std::cout << "  -i, --input <str> \t \t Input file  name containing reads (in fasta/fastq format; can be compressed)."
               " A list of files containing file names in each line can be passed with @ prefix.\n";
  std::cout << "  -k, --kmer-len <int> \t \t The length of the kmer (must be positive).\n";
  std::cout << " Optional args:\n";
  std::cout << "  -o, --output-prefix <str> \t Prefix of the output filename. [DEAFULT: SUK_k<x> where x is the value of k] \n ";
  std::cout << " -t, --threads <int> \t \t Number of threads. [DEFAULT: 1] \n";
  std::cout << "  -d, --dump-txt \t \t Dumps (in addition to the bit-vector) the unique kmers into textfile (<prefix>.txt where prefix is the value of output-prefix)."
               " [DEFAULT: False].\n";
  std::cout << "  -e, --exclude-hp \t \t Excludes the kmers having homopolymer terminals (see README.md)."
               " [DEFAULT: False].\n";
  std::cout << "  -c, --coverage <int> \t\t Expected coverage of the reads. KMC counter will be set to 4 x expected coverage."
               " [DEFAULT: 50].\n";
  std::cout << "  -m, --kmc-memory <int> \t Max memory usage of KMC (in GB)."
               " [DEFAULT: 12].\n";
  
  std::cout << "  -h, --help \t \t \t Prints the usage.\n";
}

using InputFlags = struct SInputFlags{
    std::vector<std::string> read_filenames;
    std::string output_filename;
    UINT32 threads;
    UINT16 k;
    UINT16 expected_coverage;
    UINT16 kmc_memory;
    bool exclude_hp;
    bool dump_txt; 
  };


static struct option long_options[] = {
    {"input", required_argument, NULL, 'i'},
    {"kmer-len", required_argument, NULL, 'k'},
    {"output-prefix", required_argument, NULL, 'o'},
    {"threads", required_argument, NULL, 't'},
    {"dump-txt", no_argument, NULL, 'd'},
    {"exclude-hp", no_argument, NULL, 'e'},
    {"coverage", required_argument, NULL, 'c'},
    {"kmc-memory", required_argument, NULL, 'm'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};



/** Decode the input flags
   */
void decodeFlags(int argc, char *argv[], InputFlags &flags)
{
  int args = 0;
  int opt;
  int given_k = 0;
  std::string infile;

  flags.output_filename = "SUK";
  flags.threads = 1;
  flags.expected_coverage = 50;
  flags.kmc_memory = 12;
  flags.exclude_hp = false;
  flags.dump_txt = false;

  bool is_i = false;
  bool is_k = false;

  /* initialisation */
  while ((opt = getopt_long(argc, argv, "i:k:o:t:dec:m:h", long_options,
                            nullptr)) != -1)
  {
    switch (opt)
    {
    case 'i':
      infile = optarg;
      if(infile[0]=='@') {
          std::string name;
          std::ifstream sr_list(infile.substr(1,infile.size()).c_str());
          while(std::getline(sr_list,name))
            if (name.size()>0) {
              flags.read_filenames.push_back(name);
            }
      }
      else {
          flags.read_filenames.push_back(infile);
      }
      is_i = true;
      args++;
      break;
    
    case 'k':    
      given_k = atoi(optarg);
      if (given_k > 0) {
          is_k = true;
          flags.k = (UINT16) given_k;
      }      
      args++;
      break;

    case 'o':
      flags.output_filename = std::string(optarg);
      args++;
      break;

    case 't':
      //flags.threads = std::max((UINT32)atoi(optarg)-1,(UINT32)1);
      flags.threads = std::max((UINT32)atoi(optarg),(UINT32)1);
      args++;
      break;
    case 'd':
      flags.dump_txt=true;
      break;
    case 'e':
      flags.exclude_hp=true;
      break;
    case 'c':
      flags.expected_coverage = std::max((UINT32)atoi(optarg), (UINT32)1);
      args++;
      break;
    case 'm':
      flags.kmc_memory = std::max((UINT32)atoi(optarg), (UINT32)1);
      args++;
      break;
    default:
      usage();
      exit(0);
    }
  }
  
  if (is_i && is_k )
  {
    flags.output_filename += ("_k" + std::to_string(flags.k)); 
  }
  else
  {
    fprintf(stderr, "[SUK::] Error: Invalid command: Too few arguments!\n");
    usage();
    exit(1);
  }
}

} //namespace

int main(int argc, char **argv) {

  /* Decode arguments */
  suk::InputFlags flags;
  suk::decodeFlags(argc, argv, flags);

  suk::SolidKmers sk(flags.k);
  sk.initialise(flags.read_filenames,flags.threads,flags.kmc_memory,flags.expected_coverage,flags.exclude_hp,"tmp/");
  sk.store(flags.output_filename+".bv");

  if (flags.dump_txt) {
    sk.dump_txt(flags.output_filename+".txt");
  }

  

  return 0;

}
