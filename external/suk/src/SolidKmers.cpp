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
/** Defines the class SolidKmers.
 * It represents the collection of solid kmers as sdsl bit-vector.
 */

#include <cstdlib>
#include <cassert>
#include <fstream>
#include <cstdint>
#include <stdio.h>
#include <iterator>

#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include "slog/Monitor.hpp"
#include "suk/SolidKmers.hpp"
#include "kmc_file.h"
#include "kmer_api.h"
#include <zlib.h>

namespace suk
{

SolidKmers::SolidKmers(const UINT k) : 
_k(k), _num_Solid_kmers(0),_bv(((1ULL<<(2*k))),0){}


bool SolidKmers::load(const std::string infile) {
    if (sdsl::load_from_file(_bv, infile)) {
        //_num_Solid_kmers = sdsl::sd_vector<>::rank_1_type(&_bv)(_bv.size());
        _num_Solid_kmers = sdsl::bit_vector::rank_1_type(&_bv)(_bv.size());
        assert(_bv.size() == ((1ULL<<(2*_k))));
        return true;
    }
    else {
        fprintf(stderr, "[SolidKmers] Error: Bit vector file could not be loaded: %s.\n",infile.c_str());
    }
    return false;
}

bool SolidKmers::initialise(const std::vector<std::string> & filenames, const UINT32 threads, const UINT32 max_memory, const UINT32 coverage, const bool exclude_hp, const std::string tmp_directory) {
    if(filenames.size() == 0) {
        fprintf(stderr, "[SolidKmers] Error: Empty filenames.\n");
        return false;
    }
    
    bool result = true;
    slog::Monitor monitor;
    monitor.start();
    // create tmp directory for KMC and KMC results
    struct stat st = {0};
    if (stat(tmp_directory.c_str(), &st) == -1) {
        mkdir(tmp_directory.c_str(), 0777);
    }
    std::string tmp_kmc_directory = tmp_directory + "/kmc_tmp";
    if (stat(tmp_kmc_directory.c_str(), &st) == -1) {
        mkdir(tmp_kmc_directory.c_str(), 0777);
    }
    std::string kmc_output_path = tmp_directory + "/kmc_result.res";
    // create the temp file for filenames
    std::string kmc_inputs_path = tmp_directory + "/inputs.txt";
    std::ofstream ofs(kmc_inputs_path.c_str());
    for(auto & i : filenames) {
        ofs << i << "\n";
    }
    ofs.close();
    kmc_inputs_path = "@" + kmc_inputs_path;
    
    std::string file_type = "-fq";
    gzFile fp = gzopen(filenames[0].c_str(), "r");
    char read[2];
    gzgets(fp, read, 2);
    if(read[0] == '>') {
        file_type = "-fm";
    } else if(read[0] != '@'){
        fprintf(stderr, "[SUK::] Cannot identify filetype. Assuming FASTQ for KMC input.\n");
    }
    
    // run KMC (unix only for fork(), assuming KMC is installed on machine)
    int c_pid = fork();
    if(c_pid == 0) {
        char karg[10];
        char marg[10];
        char targ[10];
        char csarg[10];
        char cxarg[10];
        char ciarg[10];
        sprintf(karg, "-k%d", _k);
        sprintf(marg, "-m%d", max_memory);
        sprintf(targ, "-t%d", threads);
        sprintf(csarg, "-cs%d", coverage*4);
        sprintf(cxarg, "-cx%d", coverage*4);
        sprintf(ciarg, "-ci%d", 2);
        //std::cout << "KMC args:" << karg << " " << marg << " " << targ << " " << csarg << " " << file_type.c_str() << " " << filename.c_str() << " " << kmc_output_path.c_str() << " " << tmp_kmc_directory.c_str() << std::endl;
        execlp("kmc", "kmc", karg, marg, targ, csarg, cxarg, ciarg, file_type.c_str(), kmc_inputs_path.c_str(), kmc_output_path.c_str(), tmp_kmc_directory.c_str(), (char *)NULL);
    }
    int status;
    wait(&status);
    if(!WIFEXITED(status)) {
        fprintf(stderr, "[SolidKmers] Error: Failed running KMC.\n");
        return false;
    }
    monitor.stop("[SUK:KMC]: Running KMC done. ");
    
    monitor.start();
    const UINT cHIST_FREQ = (coverage*4U);
    CKMCFile kmc_database;
    if(!kmc_database.OpenForListing(kmc_output_path)) {
        fprintf(stderr, "Failed to open KMC database.\n");
        return false;
    }
    
    /* Construct histogram */
    // hist0 is total kmers; hist[1] is total distinct kmers (F0); from 2, hist[i] represents number of kmers with freq=i-1
    CKMCFileInfo kmc_info;
    kmc_database.Info(kmc_info);
    
    //vector of pair kmer, counter objects
    CKmerAPI * kmer_object = new CKmerAPI(kmc_info.kmer_length);
    
    std::vector<size_t> histArray(cHIST_FREQ + 1);
    
    //UINT64 counter;
    uint64 counter;
    while (kmc_database.ReadNextKmer(*kmer_object, counter)) {
        if(counter <= cHIST_FREQ) histArray[counter]++;
    }
    monitor.stop("[SUK:KMC]: Kmers Histogram done. ");     
    
    /* Find cut offs */
    //monitor.start();
    CutOffs coffs = find_cutoffs(histArray);
    //monitor.stop("[SUK]: Finding cutoffs. ");
    
    // set bit-vector
    //monitor.start();
    UINT64 shift = 2 * (_k - 1);
    UINT64 mask = (1ULL<<2*_k) - 1;
    if(!kmc_database.RestartListing()) {
        fprintf(stderr, "Failed to restart KMC listing.\n");
        return false;
    }
    while (kmc_database.ReadNextKmer(*kmer_object, counter)) {
        if(counter >= coffs.lower && counter <= coffs.upper) {
            bool valid = true;
            if (exclude_hp && (kmer_object->get_num_symbol(0)==kmer_object->get_num_symbol(1) || kmer_object->get_num_symbol(_k-1)==kmer_object->get_num_symbol(_k-2))) { // HP at terminals; 
                valid = false;
            }
            if (valid) {
                UINT64 fwd_kmer = 0;
                UINT64 rc_kmer = 0;
                for (UINT c = 0; c < _k; c++) {
                    BYTE b = kmer_object->get_num_symbol(c);//cNt4Table[BYTE(c)];
                    if (b > 3) {
                        fprintf(stderr, "[SolidKmers] Error: Wrong base (Can not pack in 2 bits): Base %c in the kmer %s is not A, C, G, or T !\n",c,kmer_object->to_string().c_str());
                        result = false;
                        break;
                    }
                    fwd_kmer = ((fwd_kmer << 2ULL) | b) & mask;
                    rc_kmer = (rc_kmer >> 2ULL) | (3ULL^b) << shift;
                }
                _bv[fwd_kmer] = 1;
                _bv[rc_kmer] = 1;
                ++_num_Solid_kmers;
            }
        }
    }
    //monitor.stop("[SUK]: Filling bitvectors. ");

    //monitor.start();
    kmc_inputs_path = kmc_inputs_path.substr(1);
    std::string kmc_output_path_pre = kmc_output_path + ".kmc_pre";
    std::string kmc_output_path_suf = kmc_output_path + ".kmc_suf";
    unlink(kmc_inputs_path.c_str());
    //unlink(kmc_output_path_pre.c_str());
    //unlink(kmc_output_path_suf.c_str());
    rmdir(tmp_kmc_directory.c_str());
    //monitor.stop("[SUK]: Clearing files. ");

    
    fprintf(stdout, "[SolidKmers] Info: Number of solid kmers found: %lu\n",sdsl::bit_vector::rank_1_type(&_bv)(_bv.size())); 
    monitor.total("[SolidKmers]: Overall. ");
    return result;
}


bool SolidKmers::dump_txt(const std::string ofile) {
    std::ofstream ofs(ofile.c_str());
    if (!ofs.is_open()) {
        fprintf(stderr, "[SolidKmers] Error: File open error: Output file %s could not be opened !\n",ofile.c_str());
        return false;
    }
    UINT64 mask = (1ULL<<2*_k) - 1;
    UINT64 bit_mask = UINT64(0x03);
    for (UINT64 i=0; i<_bv.size(); ++i) {
        UINT64 kid = i & mask;
        std::string kmer;
        std::string rckmer;
        kmer.resize(_k,'N');
        rckmer.resize(_k,'N');
        for (INT j=_k-1; j >= 0; --j) {
            char c = cCode[kid & bit_mask];
            char rc = cRCCode[kid & bit_mask];
            kmer[j] = c;
            rckmer[_k-1-j] = rc;            
            kid  = kid >> 2ULL;
        }
        // check if canonical
        if (kmer.compare(rckmer) <= 0) {
            ofs << kmer << "\n";
        }
    }
    ofs.close();
    
    fprintf(stdout, "[SolidKmers] Info: Number of (canonical) Solid kmers found: %lu\n",_num_Solid_kmers); 
    
    return true;
}

bool SolidKmers::is_solid(std::string kmer) const {
    assert(kmer.size()==_k);
    UINT64 mask = (1ULL<<2*_k) - 1;
    UINT64 encoded_kmer = 0;
    for (auto c: kmer) {
        BYTE b = cNt4Table[BYTE(c)];
        if (b > 3) {
            fprintf(stderr, "[SolidKmers] Error: Wrong base (Can not pack in 2 bits): Base %c in the kmer %s is not A, C, G, or T !\n",c,kmer.c_str());
            return false;
        }
        encoded_kmer = ((encoded_kmer << 2ULL) | b) & mask;
    }
    return _bv[encoded_kmer];
}

CutOffs SolidKmers::find_cutoffs(const std::vector<size_t> & histArray) {
    
    /* find the cut-offs */
    // We are ignoring the last freq in hitogram as it clubs together count of kmers with freq higher than that.
    // hist0 and hit1 should be 0 (because minimum frq is set to 2 for kmc)
    CutOffs coffs;
    const int len = histArray.size()-1;
    // Ignore errornous kmers (Initial peak)
    int ind = 2;
    while (ind<len && histArray[ind] > histArray[ind+1]) {
        ++ind;
    }
    int err_th = (ind > 100)? (2):(ind); 
    coffs.err = UINT(err_th);

    // Find mean coverage value by finding global maxima)
    UINT global_maxima_val = 0;
    for (ind = err_th+1; ind<len; ++ind) {
        if (histArray[ind] > global_maxima_val) { 
            global_maxima_val = histArray[ind];
            coffs.mean = UINT(ind);
        }
    }

    // Find lower cut-off by scanning on the left of maxima (until err-th) for a freq after which most values (of nxt 5) are higher
    const int cMax_lookup = 5;
    int bind = coffs.mean-1; // ind pvs to mean
    int eind = err_th;
    // initialise
    coffs.lower = UINT(eind);
    UINT count_ge = 0;
    UINT count_lower = 0;
    for(ind = bind; ind >= eind; --ind) {
        count_ge = 0; count_lower = 0;
        for(int ind2 = ind - 1; ind2 >= (ind - cMax_lookup) && ind2 >= eind; --ind2) {
            if(histArray[ind2] < histArray[ind]) {
                ++count_lower;
            }
            else {++count_ge;}
        }
        if(count_ge >= count_lower) {
            coffs.lower = UINT(ind);
            break;
        }
    }

    // Find upper cut-off by scanning on the right of maxima for a freq after which most values (of nxt 5) are higher
    bind = coffs.mean+1; // ind nxt to mean
    eind = std::min(bind + 2*(coffs.mean-coffs.lower),UINT(len));
    // initialise
    coffs.upper = UINT(eind);
    // See if this plan works
    bool plan_a = false;
    for(ind = bind; ind < eind; ++ind) {
        count_lower = 0; count_ge = 0;
        for(UINT ind2 = ind + 1; ind2 <= ind + cMax_lookup && ind2 < len; ++ind2) {
            if(histArray[ind2] < histArray[ind]) {
                ++count_lower;                
            }
            else {++count_ge;}
        }
        if(count_ge >= count_lower) {
            coffs.upper = UINT(ind);
            plan_a = true;
            break;
        }
    }

    // Plan B: Vector of delta_avg for cases where upper threshold could not be set by this approach;
    // Stores avg delta % (delta is difference between current count and the counts in a window of next 5 which are lower than current)
    if (!plan_a){
        std::vector<UINT> delta_avg(eind,0); 
        for(ind = bind; ind < eind; ++ind) {
            UINT delta_sum=0; 
            count_lower = 0;
            for(int ind2 = ind + 1; ind2 <= ind + cMax_lookup && ind2 < len; ++ind2) {
                if(histArray[ind2] < histArray[ind]) {
                    ++count_lower;
                    delta_sum += (histArray[ind]-histArray[ind2]);
                }                
            }            
            delta_avg[ind] = UINT((delta_sum*100)/(count_lower*histArray[ind]));            
        }
        // Look for the first minimum of moving window(starting from current) avg of delta_avg
        float min_avg_avg_val = float(delta_avg[bind]);
        for(ind = bind; ind < eind; ++ind) {
            UINT window_len = std::min(cMax_lookup,eind-ind);
            UINT avg_delta_sum = 0;
            for(UINT ind2 = ind; ind2 < ind + window_len; ++ind2) {
                avg_delta_sum += delta_avg[ind2];
            }
            float avg_avg_val = float(avg_delta_sum)/float(window_len);
            // If lower than minimum so far, this is the new value
            if (avg_avg_val < min_avg_avg_val) {
                min_avg_avg_val = avg_avg_val;
                coffs.upper = UINT(ind);
            }
        }
    }


    fprintf(stdout, "[SolidKmers] Info: Error-threshold freq: %u, Lower-threshold freq: %u, Upper-threshold freq: %u, Mean-coverage: %u\n",coffs.err,coffs.lower,coffs.upper,coffs.mean); 

    return coffs;
    
}

} // namespace
