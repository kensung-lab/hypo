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

#pragma once

#include <cstdlib>
#include <cassert>
#include <fstream>
#include <cstdint>
#include <stdio.h>
#include <iterator>
//#include <tuple>
//#include <utility>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <stdexcept>
#include <sstream> 
#include <zlib.h> // for reading compressed .fq file
#include <htslib/kseq.h>
KSEQ_INIT(gzFile, gzread)

#include <htslib/sam.h>

namespace hypo{
#define DEBUG

  // Types
  using UINT = unsigned int;
  using INT = int;
  using INT8 = int8_t;
  using UINT8 = uint8_t;
  using BYTE = uint8_t;
  using INT16 = int16_t;
  using UINT16 = uint16_t;
  using INT32 = int32_t;
  using UINT32 = uint32_t;
  using INT64 = int64_t;
  using UINT64 = uint64_t;

  using ScoreParams = struct SScoreParams{
    INT8 sr_match_score;
    INT8 sr_misMatch_score;
    INT8 sr_gap_penalty; // must be negative
    INT8 lr_match_score;
    INT8 lr_misMatch_score;
    INT8 lr_gap_penalty; // must be negative

  };

  using InputFlags = struct SInputFlags{
    std::vector<std::string> sr_filenames;
    std::string sr_bam_filename;
    std::string lr_bam_filename;
    std::string draft_filename;
    std::string output_filename;
    ScoreParams score_params;
    UINT8 map_qual_th; // mapping_qual_th
    UINT64 norm_edit_th; // normalised edit distance threshold for long reads
    UINT32 threads;
    UINT32 processing_batch_size;
    UINT16 k;
    UINT16 cov;
    UINT16 sz_in_gb;
    UINT done_stage;
    // whether intermediate results should be used and stored to disk; 
    // if false, neither intermed files will be written nor used (even if exist)
    bool intermed; 

  };

  // Stages
  #define STAGE_BEG 0u
  #define STAGE_SK 1u
  #define STAGE_SP 2u

  // Region types (Needed by Contig and Alignment Classes)
  enum class RegionType: UINT8 {
      SWS,
      SW,
      WS,
      MWM,
      MW,
      WM,
      SWM,
      MWS,
      OTHER,
      LONG,
      SR,
      MSR // Minimser is same as SR (no polishing)
  };
  

  // Folders and files
  #define AUX_DIR  "aux/"
  #define SKFILE "aux/solid_kmers.bvsd"
  #define STAGEFILE "aux/stage.txt"
  #define INSPECTFILEPREF "aux/inspect_"
  #define BEDFILE "aux/regions.bed"

  // VArious settings
  #define SR_COV_TH 5u
  #define SR_SUPP_FRAC 0.4
  //#define SR_SUPP_FRAC 0.8
  // Minimiser based window cutting thresholds/constants
  #define MINIMIZER_K 10u // Change Poly-base macros, if this k gets changed.
  #define MINIMIZER_W 10u
  #define MINIMISER_COV_TH 5u
  #define MINIMISER_SUPP_FRAC 0.8
  #define IDEAL_WIND_SIZE 100u
  #define WIND_SIZE_SEARCH_TH 80u //This should always be <= ideal_wind_size
  #define POLYA 0x000000u
  #define POLYC 0x055555u
  #define POLYG 0x0aaaaau
  #define POLYT 0x0fffffu

  //#define IDEAL_WIND_SIZE 140u
  #define IDEAL_LWIND_SIZE 500u // Ideal long window size
  //#define WIND_SIZE_SEARCH_TH 120u //This should always be <= ideal_wind_size
  #define MIN_SHORT_NUM 3 // Minimum num of short-reads arms for a window to be declared short  
  #define MIN_INTERNAL_NUM1 20u // Minimum num of internal arms for a short-arms window to discard pref/suff arms
  #define MIN_INTERNAL_NUM2 5u // Minimum num of internal arms for a spcl (SW,WS,SWS) short-arms window to discard pref/suff arms
  #define MIN_INTERNAL_NUM3 10u // Minimum num of internal arms for a long-arms window to discard pref/suff arms
  #define MIN_CONTRIB 10u // Minimum num of arms for a short-arms window to be considered for discarding pref/suff arms
  #define MIN_INTERNAL_CONTRIB 0.4 // Minimum fraction of internal arms (wrt total arms) for a short-arms window to discard pref/suff arms
  #define SHORT_ARM_COEF 10 // Short arm len should be >= window_len/coeffcient



  const std::vector<char> cCode = {'A','C','G','T','N'};
// A: 0, C: 1, G:2, T: 3; everything else is 4.
// Unpacking will map everything other ACGT to N.
const BYTE cNt4Table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
  
  enum Stage{
    SK, // solid kmer
    SP // solid pos

      };
} // end namespace


