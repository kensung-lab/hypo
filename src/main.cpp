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

#include <sys/time.h>
#include <getopt.h>
#include <cctype>
#include <math.h> 
#include <algorithm>
#include <sys/stat.h>

#include "globalDefs.hpp"
#include "Hypo.hpp"


/** Module containing main() method for reading and processing arguments.
 */

namespace hypo{
void usage (void);
void decodeFlags(int argc, char* argv [], InputFlags& flags);
std::vector<std::string> split(const std::string &str, char delimiter);
void print_byte_to_strings(int nb);
UINT get_kmer_len(const std::string& arg);
UINT16 get_expected_file_sz(const std::string& given_size, UINT16 cov);
void set_kind(const std::string& kind);

static struct option long_options[] = {
    {"reads-short", required_argument, NULL, 'r'},
    {"draft", required_argument, NULL, 'd'},
    {"size-ref", required_argument, NULL, 's'},
    {"coverage-short", required_argument, NULL, 'c'},
    {"bam-sr", required_argument, NULL, 'b'},
    {"bam-lr", required_argument, NULL, 'B'},
    {"output", required_argument, NULL, 'o'},
    {"threads", required_argument, NULL, 't'},
    {"processing-size", required_argument, NULL, 'p'},
    {"kind-sr", required_argument, NULL, 'k'},
    {"match-sr", required_argument, NULL, 'm'},
    {"mismatch-sr", required_argument, NULL, 'x'},
    {"gap-sr", required_argument, NULL, 'g'},
    {"match-lr", required_argument, NULL, 'M'},
    {"mismatch-lr", required_argument, NULL, 'X'},
    {"gap-lr", required_argument, NULL, 'G'},
    {"qual-map-th", required_argument, NULL, 'q'},
    {"ned-th", required_argument, NULL, 'n'},
    {"intermed", no_argument, NULL, 'i'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};

inline bool file_exists (const std::string& name) {
  struct stat st;   
  return (stat (name.c_str(), &st) == 0); 
}

inline bool create_dir (const std::string& name) {
  struct stat st; 
  int status = 0;  
  if (stat(name.c_str(), &st) == -1) {
    status = mkdir(name.c_str(), 0777);
  }
  return (status==0);
}

/** Define various settings as globals
  */
const SRSettings Sr_settings = {5u,0.4};
const MinimizerSettings Minimizer_settings = {10u,10u,5u,0.8,0x000000u,0x055555u,0x0aaaaau,0x0fffffu};
WindowSettings Window_settings = {100u,500u,80u};
const ArmsSettings Arms_settings = {3u,20u,5u,10u,10u,0.4,10u};


/** Decode the input flags
   */
void decodeFlags(int argc, char *argv[], InputFlags &flags)
{
  int args = 0;
  int opt;
  std::string infile;

  std::string kind="sr";
  flags.lr_bam_filename = "";
  flags.score_params.sr_match_score = 5;
  flags.score_params.sr_misMatch_score = -4;
  flags.score_params.sr_gap_penalty = -8;
  flags.score_params.lr_match_score = 3;
  flags.score_params.lr_misMatch_score = -5;
  flags.score_params.lr_gap_penalty = -4;
  flags.map_qual_th = 2;
  flags.norm_edit_th = 20;
  flags.threads = 1;
  flags.processing_batch_size = 0;
  flags.output_filename = "";
  flags.intermed = false;
  flags.sz_in_gb = 12;

  bool is_sr = false;
  bool is_draft = false;
  bool is_size = false;
  bool is_cov = false;
  bool is_bamsr = false;
  std::string given_sz;
  std::string cmd = "hypo ";
  std::string err_string = "";
  /* initialisation */
  while ((opt = getopt_long(argc, argv, "r:d:s:c:b:B:o:t:p:k:m:x:g:M:X:G:q:n:ih", long_options,
                            nullptr)) != -1)
  {
    switch (opt)
    {
    case 'r':
      infile = optarg;
      cmd += (" -r " + std::string(optarg));
      if(infile[0]=='@') {
          std::string name;
          std::ifstream sr_list(infile.substr(1,infile.size()).c_str());
          if (!sr_list.good()) {
              fprintf(stderr, "[Hypo::utils] Error: File Error: Could not open the file %s!\n",infile.substr(1,infile.size()).c_str());
              exit(1); 
          }
          while(std::getline(sr_list,name)) {
            if (name.size()>0) {
              flags.sr_filenames.push_back(name);
              if (!file_exists(name)) {
                fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",name.c_str());
                exit(1);
              }
            }
        }
      }
      else {
          flags.sr_filenames.push_back(infile);
          if (!file_exists(infile)) {
            fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",infile.c_str());
            exit(1);
          }
      }
      is_sr = true;
      args++;
      break;
    
    case 'd':
      flags.draft_filename = std::string(optarg);
      if (!file_exists(flags.draft_filename)) {
        fprintf(stderr, "[Hypo::utils] Error: File Error: Draft file does not exist %s!\n",flags.draft_filename.c_str());
        exit(1);
      }
      cmd += (" -d " + std::string(optarg));
      args++;
      is_draft = true;
      break;

    case 's':    
      flags.k = std::max(2U,get_kmer_len(std::string(optarg)));
      cmd += (" -s " + std::string(optarg));
      given_sz = std::string(optarg);
      is_size = true;
      args++;
      break;

    case 'c': 
      if (atoi(optarg) <= 0 ) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Coverage should be positive %d!\n",atoi(optarg));
        exit(1);
      }   
      flags.cov = atoi(optarg);
      cmd += (" -c " + std::string(optarg));
      is_cov = true;
      args++;
      break;

    case 'b':
      flags.sr_bam_filename = std::string(optarg);
      if (!file_exists(flags.sr_bam_filename)) {
        fprintf(stderr, "[Hypo::utils] Error: File Error: Short reads BAM file does not exist %s!\n",flags.sr_bam_filename.c_str());
        exit(1);
      }
      cmd += (" -b " + std::string(optarg));
      is_bamsr = true;
      args++;
      break;

    case 'B':
      flags.lr_bam_filename = std::string(optarg);
      if (!file_exists(flags.lr_bam_filename)) {
        fprintf(stderr, "[Hypo::utils] Error: File Error: Long reads BAM file does not exist %s!\n",flags.lr_bam_filename.c_str());
        exit(1);
      }
      cmd += (" -B " + std::string(optarg));
      args++;
      break;

    case 'o':
      flags.output_filename = std::string(optarg);
      cmd += (" -o " + std::string(optarg));
      args++;
      break;

    case 'k':
      kind = std::string(optarg);
      cmd += (" -k " + kind);
      args++;
      break;

    case 'm':
      flags.score_params.sr_match_score = atoi(optarg);
      cmd += (" -m " + std::string(optarg));
      args++;
      break;
    case 'x':
      flags.score_params.sr_misMatch_score = atoi(optarg);
      cmd += (" -x " + std::string(optarg));
      args++;
      break;
    case 'g':
      if (atoi(optarg) >= 0 ) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Gap penalty (g) must be negative %d!\n",atoi(optarg));
        exit(1);
      } 
      flags.score_params.sr_gap_penalty = atoi(optarg);
      cmd += (" -g " + std::string(optarg));
      args++;
      break;
    case 'M':
      flags.score_params.lr_match_score = atoi(optarg);
      cmd += (" -M " + std::string(optarg));
      args++;
      break;
    case 'X':
      flags.score_params.lr_misMatch_score = atoi(optarg);
      cmd += (" -X " + std::string(optarg));
      args++;
      break;
    case 'G':
      if (atoi(optarg) >= 0 ) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Gap penalty (G) must be negative %d!\n",atoi(optarg));
        exit(1);
      }
      flags.score_params.lr_gap_penalty = atoi(optarg);
      cmd += (" -G " + std::string(optarg));
      args++;
      break;
    case 'q':
      if (atoi(optarg) < 0 ) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Mapping quality threshold (q) must NOT be negative %d!\n",atoi(optarg));
        exit(1);
      }
      flags.map_qual_th = atoi(optarg);
      cmd += (" -q " + std::string(optarg));
      args++;
      break;
    case 'n':
      if (atoi(optarg) < 0) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Normalised Edit Distance Threshold (n) must NOT be negative %d!\n",atoi(optarg));
        exit(1);
      }
      flags.norm_edit_th = atoi(optarg);
      cmd += (" -n " + std::string(optarg));
      args++;
      break;
    case 't':
      if (atoi(optarg) <= 0) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Number of threads (t) must be positive %d!\n",atoi(optarg));
        exit(1);
      }
      //flags.threads = std::max((UINT32)atoi(optarg)-1,(UINT32)1);
      flags.threads = std::max((UINT32)atoi(optarg),(UINT32)1);
      cmd += (" -t " + std::string(optarg));
      args++;
      break;
    case 'p':
      if (atoi(optarg) <= 0) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Processing-size, i.e. number of contigs processed in a batch, (p) must NOT be negative %d!\n",atoi(optarg));
        exit(1);
      }
      flags.processing_batch_size = std::max((UINT32)atoi(optarg),(UINT32)0);
      cmd += (" -p " + std::string(optarg));
      args++;
      break;
    case 'i':
      flags.intermed=true;
      cmd += (" -i ");
      break;
    default:
      usage();
      exit(0);
    }
  }
  // TODO: Check if the int args conform to the assumpions
  
  if (is_sr && is_draft && is_size && is_bamsr && is_cov)
  {
    // Set window settings
    void set_kind(const std::string& kind);
    // Set expected short reads file size
    flags.sz_in_gb = get_expected_file_sz(given_sz,flags.cov);

    // Set output name
    if (flags.output_filename == "") {
      size_t ind = flags.draft_filename.find_last_of("(/\\");
      std::string dfullname (flags.draft_filename,ind + 1);
      ind = dfullname.find_last_of(".");
      std::string dname(dfullname,0,ind);
      flags.output_filename = "hypo_" + dname + ".fasta";
    }
    fprintf(stdout, "Given Command: %s.\n",cmd.c_str()); 
    // Set stage to start from  
    create_dir(AUX_DIR);
    if (flags.intermed) {  
      if (file_exists(STAGEFILE)) {
        std::ifstream ifs(STAGEFILE);
        if (!ifs.is_open()) {
          fprintf(stderr, "[Hypo::Utils] Error: File open error: Stage File (%s) exists but could not be opened!\n",STAGEFILE);
          exit(1);
        }
        std::string dummy1,dummy2,dummy3;
        UINT stage_num=STAGE_BEG;
        while (ifs >> dummy1 >> dummy2 >> dummy3 >> stage_num){}
        flags.done_stage = stage_num;
        std::cout << "Stagenum found is " << stage_num <<std::endl;
        // TODO: Check whether all files required for later stages are present
      }
      else {
        flags.done_stage = STAGE_BEG;
      }
      fprintf(stdout, "[Hypo::Utils] Info: Intermediate Files will be stored.\n"); 
    }
    else {
      flags.done_stage = STAGE_BEG;
      fprintf(stdout, "[Hypo::Utils] Info: Intermediate Files will NOT be stored.\n"); 
    }
    fprintf(stdout, "[Hypo::Utils] Info: Beginning from stage: %u\n",flags.done_stage); 
  }
  else
  {
    fprintf(stderr, "[Hypo::] Error: Invalid command: Too few arguments!\n");
    usage();
    exit(1);
  }
}

/*
   * Usage of the tool
   */
void usage(void)
{
  std::cout << "\n Usage: hypo <args>\n\n";
  std::cout << " ** Mandatory args:\n";
  std::cout << "\t-r, --reads-short <str>\n"
            << "\tInput file name containing reads (in fasta/fastq format; can be compressed). "
            << "A list of files containing file names in each line can be passed with @ prefix.\n\n";
  std::cout << "\t-d, --draft <str>\n"
            << "\tInput file name containing the draft contigs (in fasta/fastq format; can be compressed). \n\n";
  std::cout << "\t-b, --bam-sr <str>\n"
            << "\tInput file name containing the alignments of short reads against the draft (in bam/sam format; must have CIGAR information). \n\n";
  std::cout << "\t-c, --coverage-short <int>\n"
            << "\tApproximate mean coverage of the short reads. \n\n";
  std::cout << "\t-s, --size-ref <str>\n"
            << "\tApproximate size of the genome (a number; could be followed by units k/m/g; e.g. 10m, 2.3g). \n\n\n";


  std::cout << " ** Optional args:\n";
  std::cout << "\t-B, --bam-lr <str>\n"
            << "\tInput file name containing the alignments of long reads against the draft (in bam/sam format; must have CIGAR information). \n"
            << "\t[Only Short reads polishing will be performed if this argument is not given]\n\n";
  std::cout << "\t-o, --output <str>\n"
            << "\tOutput file name. \n"
            << "\t[Default] hypo_<draft_file_name>.fasta in the working directory.\n\n ";

  std::cout << "\t-t, --threads <int>\n"
            << "\tNumber of threads. \n"
            << "\t[Default] 1.\n\n ";
  std::cout << "\t-p, --processing-size <int>\n"
            << "\tNumber of contigs to be processed in one batch. Lower value means less memory usage but slower speed. \n"
            << "\t[Default] All the contigs in the draft.\n\n ";
  std::cout << "\t-k, --kind-sr <str>\n"
            << "\tKind of the short reads. \n"
            << "\t[Valid Values] \n"
            << "\t\tsr\t(Corresponding to NGS reads like Illumina reads) \n"
            << "\t\tccs\t(Corresponding to HiFi reads like PacBio CCS reads) \n"
            << "\t[Default] sr.\n\n ";

  std::cout << "\t-m, --match-sr <int>\n"
            << "\tScore for matching bases for short reads. \n"
            << "\t[Default] 5.\n\n ";
  std::cout << "\t-x, --mismatch-sr <int>\n"
            << "\tScore for mismatching bases for short reads. \n"
            << "\t[Default] -4.\n\n ";
  std::cout << "\t-g, --gap-sr <int>\n"
            << "\tGap penalty for short reads (must be negative). \n"
            << "\t[Default] -8.\n\n ";
  std::cout << "\t-M, --match-lr <int>\n"
            << "\tScore for matching bases for long reads. \n"
            << "\t[Default] 3.\n\n ";
  std::cout << "\t-X, --mismatch-lr <int>\n"
            << "\tScore for mismatching bases for long reads. \n"
            << "\t[Default] -5.\n\n ";
  std::cout << "\t-G, --gap-lr <int>\n"
            << "\tGap penalty for long reads (must be negative). \n"
            << "\t[Default] -4.\n\n ";
  std::cout << "\t-n, --ned-th <int>\n"
            << "\tThreshold for Normalised Edit Distance of long arms allowed in a window (in %). Higher number means more arms allowed which may slow down the execution.\n"
            << "\t[Default] 20.\n\n ";
  std::cout << "\t-q, --qual-map-th <int>\n"
            << "\tThreshold for mapping quality of reads. The reads with mapping quality below this threshold will not be taken into consideration. \n"
            << "\t[Default] 2.\n\n ";
  std::cout << "\t-i, --intermed\n"
            << "\tStore or use (if already exist) the intermediate files. \n"
            << "\t[Currently, only Solid kmers are stored as an intermediate file.].\n\n ";
  std::cout << "\t-h, --help\n"
            << "\tPrint the usage. \n\n";
}


std::vector<std::string> split(const std::string &str, char delimiter) {
  std::string token;
  std::istringstream ss{str};
  std::vector<std::string> tokens;
  while (std::getline(ss, token, delimiter))
  {
    tokens.push_back(token);
  }
  return tokens;
}

void print_byte_to_strings(int nb) {
    // 0x3 for 2-bits and 0x7 for 4-bits
    BYTE base_mask = (BYTE(1) << nb)-1;
    
    // for 2-bits, every BYTE is valid
    if (nb==2) {
        for (UINT i=0; i < 256; ++i) {
            auto num = BYTE(i);
            if (i%16==0) {
              std::cout << "\n\t";
            }
            std::string token(4,'A');
            for (UINT j=4; j > 0; --j) {
                token[j-1] = cCode[(num & base_mask)];
                num = num >> nb;
            }
            if (i%4==0 && i%16!=0) {
              std::cout << " ";
            }
            std::cout << "\"" << token << "\"" << ", ";
        }
    }
    else { // for 4-bits, not every UINT8 is valid
        std::string token("NN");
        for (UINT i=0; i < 256; ++i) {
            auto num = BYTE(i);
            if (i%16==0) {
              std::cout << "\n\t";
            }
            BYTE second_qbit = num & base_mask;
            if (second_qbit < cCode.size() ) { // A,C,G,T,N
                token[1] = cCode[second_qbit];
            }
            num = num >> nb;
            BYTE first_qbit = num & base_mask;
            if (first_qbit < cCode.size()) { // A,C,G,T,N
                token[0] = cCode[first_qbit];
            }
            if (i%4==0 && i%16!=0) {
              std::cout << " ";
            }
            std::cout << "\"" << token << "\"" << ", ";
        }
    }       
}

UINT get_kmer_len(const std::string& given_size) {
  // Find minimum k such that it is odd and  4^k >= genome_size => 2**2k >= genome_size => 2k = log_2(genome_size)
  UINT power;
  size_t ind;
  auto val = std::stof(given_size,&ind);
  if (ind >= given_size.size()) { // should be absolute number
    if (floor(val)!=ceil(val)) {
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Genome-size with no units (K,M,G etc,) should be absolute number!\n");
      exit(1);
    }
    power = 0;
  }
  else {
    char unit = std::toupper(given_size[ind]);
    switch (unit)
    {
    case 'K':
      power = 10;
      break;
    case 'M':
      power = 20;
      break;
    case 'G':
      power = 30;
      break;
    case 'T':
      power = 40;
      break;
    default:
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Allowed units for Genome-size are K (10^3),M (10^6),G (10^9),T (10^12)!\n");
      exit(1);
    }
  }
  UINT kmer_len = (power) + ceil(log2(val));
  kmer_len = ceil(kmer_len/2); // divide by 2
  if (kmer_len%2==0) {++kmer_len;}
  fprintf(stdout, "[Hypo::Utils] Info: Value of K chosen for the given genome size (%s): %u\n",given_size.c_str(),kmer_len); 
  return kmer_len;
}

UINT16 get_expected_file_sz(const std::string& given_size, UINT16 cov) {
  // Find minimum k such that it is odd and  4^k >= genome_size => 2**2k >= genome_size => 2k = log_2(genome_size)
  UINT64 sz_in_gb=1024;
  size_t ind;
  auto val = std::stof(given_size,&ind);
  val = 2 * cov * val;
  if (ind >= given_size.size()) { // should be absolute number
    if (floor(val)!=ceil(val)) {
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Genome-size with no units (K,M,G etc,) should be absolute number!\n");
      exit(1);
    }
    else {
      sz_in_gb = val/1000000000;
    }
  }
  else {
    char unit = std::toupper(given_size[ind]);
    switch (unit)
    {
    case 'K':
      sz_in_gb = val/1000000;
      break;
    case 'M':
      sz_in_gb = val/1000;
      break;
    case 'G':
      sz_in_gb = val;
      break;
    case 'T':
      sz_in_gb = 1024;
      break;
    default:
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Allowed units for Genome-size are K (10^3),M (10^6),G (10^9),T (10^12)!\n");
      exit(1);
    }
  }
  if (sz_in_gb < 12) {sz_in_gb=12;}
  else if (sz_in_gb > 1024) {sz_in_gb=1024;}
  fprintf(stdout, "[Hypo::Utils] Info: File size expected for the given genome size (%s) and cov (%u): %luG\n",given_size.c_str(),cov,sz_in_gb); 
  return UINT16(sz_in_gb);
}

void set_kind(const std::string& kind) {
  if (kind == "sr") {
    Window_settings.ideal_swind_size = 100;
    Window_settings.wind_size_search_th = 80;
  }
  else if (kind == "ccs") {
    Window_settings.ideal_swind_size = 500;
    Window_settings.wind_size_search_th = 400;
  }
  else {
    fprintf(stderr, "[Hypo::Utils] Error: Wrong value for kind-sr: Allowed values are <sr>, <ccs>!\n");
    exit(1);
  }
}


} // end namespace

int main(int argc, char **argv) {

  /* Decode arguments */
  hypo::InputFlags flags;
  hypo::decodeFlags(argc, argv, flags);

  hypo::Hypo hypo(flags);
  hypo.polish();
  return 0;
}
