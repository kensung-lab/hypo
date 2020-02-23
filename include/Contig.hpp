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
/** Class Contig.
 * It represents a contig and contains the functionality for different processing that will be carried out on a contig.
 */
#pragma once
#include <memory>
#include <mutex>
#include <unordered_map>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include "suk/SolidKmers.hpp"

#include "globalDefs.hpp"
#include "PackedSeq.hpp"
#include "Alignment.hpp"
#include "Window.hpp"

namespace hypo
{
#define MAX_U16_LIMIT 0xffffu
using KmerInfo = struct _kmerinfo {
    _kmerinfo(UINT64 kmer):kid(kmer),coverage(0),support(0){}
    UINT64 kid;
    UINT16 coverage;
    UINT16 support;
    std::mutex k_mutex;
};

using MWMinimiserInfo = struct Sminimiserinfo {
    std::vector<UINT32> minimisers; // minimser (2 bit based); The MSB represents validity (0 for valid, 1 for invalid)
    std::vector<UINT32> rel_pos; // Relative position of each minimiser wrt pvs one' first is wrt the start of the window
    std::vector<UINT16> support;
    std::vector<UINT16> coverage;
    std::mutex m_mutex;
};

using DivisionPoint = struct Sdivisionpoint {
    UINT32 start_pos; // position of the first base of the window
    UINT32 minimiser;
};



class Contig {
public:
    Contig(const UINT32 id, const std::string& name, const std::string seq);
    Contig(const UINT32 id, const kseq_t * ks);
    ~Contig() = default;
    // Fills _solid_pos and _kmerinfo
    void find_solid_pos(const std::unique_ptr<suk::SolidKmers>& pSK);

    // This will prepare the contig for minimer based window-divsion; called after kmer-support has been updated.
    void prepare_for_division(const UINT k);

    // This will divide the contig into regions; called after minimser-support has been updated.
    void divide_into_regions();
    
    /** This will first fill arms into windows.
     * Then curate the short windows; resetting thosewith insufficient short-arms; 
     * Called after preparing arms.
     * Will destroy the alignments after use
     * */
    void fill_short_windows(std::vector<std::unique_ptr<Alignment>>& alignments);
    
    // This will set the ds associated with long-windows; Called before long-arm filling.
    void prepare_long_windows();


    /** This will first fill arms into long windows
     * Then destroy ds associated with long-windows; Called after long-arm preparing.
     * Alignments will be destroyed
     * */
    inline void fill_long_windows(std::vector<std::unique_ptr<Alignment>>& alignments) {
        /* Fill Short Windows */
        for (UINT i=0; i < alignments.size();++i) {
            alignments[i]->add_arms(*this);
            alignments[i].reset();
        }
        // Get rid of prefix-suffix in long windows, if possible
        auto num_reg = _reg_type.size()-1; //excluding the dummy
        for (UINT32 i=0; i < num_reg; ++i) {
            if (_reg_type[i] == RegionType::LONG) { // a long window
                if (_pwindows[i]->get_num_internal()>Arms_settings.min_internal_num3) {
                    _pwindows[i]->clear_pre_suf();
                }                
            }
        }
        sdsl::util::clear(_pseudo_reg_pos);
        sdsl::util::clear(_pseudo_Rreg_pos);
        sdsl::util::clear(_pseudo_Sreg_pos);
        _pseudo_reg_type.clear();
        _pseudo_reg_type.shrink_to_fit();
        _true_reg_id.clear();
        _true_reg_id.shrink_to_fit();
    }

    inline UINT64 get_num_regions() const {return _reg_type.size()-1;}
    inline UINT64 get_num_sr() const {return _numSR;}
    inline UINT64 get_len_sr() const {return _lenSR;}
    inline bool is_valid_window(UINT32 ind) {assert(ind<_reg_type.size()); return (_pwindows[ind]!=nullptr);}
    inline void generate_consensus(const UINT64 ind, const UINT32 th) {_pwindows[ind]->generate_consensus(th);} 

    void generate_inspect_file(std::ofstream& bedfile);

    static void set_no_long_reads() {_no_long_reads = true;}

    friend std::ostream &operator<<(std::ostream &, const Contig &);
    friend class Alignment;

    

private:
    const UINT32 _id;
    const std::string _name;
    const UINT32 _len;
    PackedSeq<4> _pseq;

    // Created by find_solid_pos; Will be used for finding SR; Will be discarded after that
    sdsl::bit_vector _solid_pos;
    sdsl::bit_vector::rank_1_type _Rsolid_pos;
    sdsl::bit_vector::select_1_type _Ssolid_pos;
    std::vector<std::unique_ptr<KmerInfo>> _kmerinfo; 
    

    /** Store first and the last kmer-id for each SR (starting from index 1). Index 0 is dummy.
     * Every i^th SR (1-based) will have its first and the last KID stored at (2i-1) and (2i) .
     * Each window of type SWS, SW, WS, SWM, or MWS will make use of kmer in this.
     * Will be discarded after filling Short arms
     * 
     * */
    std::vector<UINT64> _anchor_kmers;

    /** DS used for dividing MegaWindows (Window between SRs) into windows
     * _reg_pos: bit-vector marking the beginning of each region (SR or MegaWindow); has a dummy at the end (= len of contig)
     *      - Later when MW will be divided into windows, more bits will be set in it (corresponding to each window and Minimiser-SR)
     * rank-select: supporting DS
     * 
     * _is_win_even: Marks whether MW are even numbered(0,2,4,..) or not (1,3,5,..)
     *   - We do not need to save type of each region as it can be deduced.
     *     * if (_is_win_even): i%2==0 => MW; otherwise => SR
     *     * if (!_is_win_even): i%2==0 => SR; otherwise => MW
     *   - Index of Hashtable for each MW can also be deduced
     *     * if (_is_win_even): i/2 is the index
     *     * Else: (i-1)/2 is the index
     * 
     * _minimserinfo: Vector containing, for each MW, its minimser info (hastable for support and corresponding mutex)
     * 
     * Except _reg_pos, rest will be destroyed after division
     * */
    sdsl::bit_vector _reg_pos;
    sdsl::bit_vector::rank_1_type _RMreg_pos;
    sdsl::bit_vector::select_1_type _SMreg_pos;
    bool _is_win_even;
    std::vector<std::unique_ptr<MWMinimiserInfo>> _minimserinfo;

    /** Together, the following 3 vectors represent region-division map of the contig
     * A region is either an SR or a window
     * 
     * _reg_info: Stores 32-bit unsigned integer showing minimser value for MSR or rank of SR (for SR)
     *      Other Region types (i.e. windows) have invalid values (set to 0);
     *      rank of SR: (i^th (1-based) in the contig SR will have i stored as info)
     *      This is done so that SW, WS , SWS, SWM, MWS know their preceding/succeeding  SR and hence the anchor-KID (through _anchor_kmers)
     *      Will be discarded after filling Short arms
     * _pwindows:
     * _pwindows will be a nullptr for SR and invalid windows 
     * (windows will become invalid when merged for Long reads)
     * */
    
    sdsl::bit_vector::rank_1_type _Rreg_pos;
    sdsl::bit_vector::select_1_type _Sreg_pos;
    std::vector<RegionType> _reg_type;
    std::vector<UINT32> _reg_info;
    std::vector<std::unique_ptr<Window>> _pwindows;
    // Total Number and the length of SR
    UINT64 _numSR;
    UINT64 _lenSR;

    // Pseudo-regions used for long reads; will be deleted after long-arm filling
    sdsl::bit_vector _pseudo_reg_pos;
    sdsl::bit_vector::rank_1_type _pseudo_Rreg_pos;
    sdsl::bit_vector::select_1_type _pseudo_Sreg_pos;
    std::vector<RegionType> _pseudo_reg_type;
    std::vector<UINT32> _true_reg_id;

    // Will be set to true if only short reads polishing; Needed while writing the results
    static bool _no_long_reads;
    
    
    inline void increment_support(const UINT32 ind) { std::lock_guard<std::mutex> guard(_kmerinfo[ind]->k_mutex); ++(_kmerinfo[ind]->support);}
    inline void increment_coverage(const UINT32 ind) { std::lock_guard<std::mutex> guard(_kmerinfo[ind]->k_mutex); ++(_kmerinfo[ind]->coverage);}
    inline void increment_minimser_support(const UINT32 mininfo_idx, const UINT32 mindex) { 
        {std::lock_guard<std::mutex> guard(_minimserinfo[mininfo_idx]->m_mutex); ++(_minimserinfo[mininfo_idx]->support[mindex]);}
    }
    inline void increment_minimser_coverage(const UINT32 mininfo_idx, const UINT32 mindex) { 
        {std::lock_guard<std::mutex> guard(_minimserinfo[mininfo_idx]->m_mutex); ++(_minimserinfo[mininfo_idx]->coverage[mindex]);}
    }
    void initialise_minimserinfo(const std::string& draft_seq, UINT32 minfoind);
    void divide(const UINT32 reg_index, const UINT32 beg, const UINT32 end, char pvs, char nxt);
    void force_divide(const UINT32 beg, const UINT32 end, char pvs, char nxt);
}; // Polisher
} // namespace hypo