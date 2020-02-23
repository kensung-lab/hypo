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
/** Defines the class Alignment.
 * It represents an alignment and contains functionality for kmer -support update; minimizer-support update, and arms finding
 */
#include <cmath>
#include "Alignment.hpp"
#include "Contig.hpp"
#include "MinimizerDeque.hpp"

namespace hypo
{
Alignment::Alignment(Contig& contig, bam1_t *hts_align): is_valid(true) {
    initialise_pos(hts_align);
    UINT32 clen = contig._len;
    if (_rb >= clen || _re > clen) {
        _qname = std::string(bam_get_qname(hts_align), hts_align->core.l_qname);
        fprintf(stdout, "[Hypo::Alignment] Error: Alignment File error: Looks like the reference in the alignment file is different from the draft. Contig (%s): Read (%s): rb (%u): re (%u): clen (%u)\n",contig._name.c_str(), _qname.c_str(), _rb, _re, clen); 
        exit(1);
    }
    copy_data(hts_align);
}

Alignment::Alignment(Contig& contig, UINT64 norm_edit_th, bam1_t *hts_align): is_valid(true) {
    initialise_pos(hts_align);
    UINT32 rlen = _re-_rb;
    UINT32 clen = contig._len;
    if (_rb >= clen || _re > clen) {
        _qname = std::string(bam_get_qname(hts_align), hts_align->core.l_qname);
        fprintf(stdout, "[Hypo::Alignment] Error: Alignment File error:: Looks like the reference in the alignment file is different from the draft. Contig (%s): Read (%s): rb (%u): re (%u): clen (%u)\n",contig._name.c_str(), _qname.c_str(), _rb, _re, clen); 
        exit(1);
    }
    //UIN32 qlen = _qae-_qab;
    //UINT32 alen = std::max(rlen,alen);
    auto nmp = bam_aux_get(hts_align,"NM");
    if (nmp) {
        INT64 edit_dist = bam_aux2i(nmp);
        auto norm_edit_dist = std::ceil(edit_dist*100/rlen);
        if (norm_edit_dist > norm_edit_th) { // edit dist above tolerance
            is_valid = false;
        }
    }
    if (is_valid) {
        copy_data(hts_align);
    }

}

void Alignment::update_solidkmers_support (const UINT k, Contig& contig) {
    /* Find kmer range */
    auto first = contig._Rsolid_pos(_rb);
    auto last = contig._Rsolid_pos(_re);
    // discard those which do not wholly fall in the alignment
    // NOTE SELECT query's index is 1-based
    for (auto i = last; i >first; --i) {
        auto pos = contig._Ssolid_pos(i);
        if (pos + k <= _re ) { // falls wholly within
            last = i;
            break;
        }
    }
    // Here first+1 th and last th pos fall in the read; index falling in the read: first and last-1; lst points to index of kmer outside (or not fully in)
    if  (last>first) { // some kmer found
        std::unordered_multimap<UINT64,UINT32> kmap;
        kmap.reserve(last-first);
        std::vector<UINT64> spos;
        spos.reserve(last-first);
        // update coverage; Add kid to map; add pos of solid kmer to spos
        for (auto i = first; i <last; ++i) {
            contig.increment_coverage(i);
            kmap.insert({contig._kmerinfo[i]->kid,i-first});
            spos.emplace_back(contig._Ssolid_pos(i+1));
        }
        // Find kmers in the read
        UINT64 kmer=0;
        UINT kmer_len=0;
        const UINT64 kmask = (1ULL<<2*k) - 1;
        size_t num_qbases = _apseq.get_seq_size();
        UINT32 num_cbases = _re-_rb;
        INT64 pvs_supp_kpos = -1; 
        UINT32 pvs_supp_r_bind = 0;
        for (size_t r_ind=0; r_ind < num_qbases; ++r_ind) { //r_ind is index of the end of kmer
            BYTE b = _apseq.enc_base_at(r_ind);
            //Always ACGT for reads
            kmer = (kmer << 2 | b) & kmask;
            if (kmer_len<k) { ++kmer_len;}   
            // Update support if conditions meet           
            if (kmer_len==k) {
                UINT32 r_bind = r_ind+1-k;
                auto kit = kmap.equal_range(kmer);
                for  (auto itr=kit.first; itr!=kit.second;++itr) {
                    UINT32 c_ind = itr->second;
                    INT64 c_dist = spos[c_ind] - _rb;
                    UINT32 srange_left = ((c_dist>k)?(c_dist-k):(0));
                    UINT32 srange_right = std::min(INT64(num_cbases),c_dist+(k));
                    // WARNING: The same read-kmer can provide support to many in the contig
                    // (if it is not unique and repeat close to each other). But such cases  are rare. So ignore.
                    if (r_bind >= srange_left && r_bind <= srange_right) { // supports
                        bool should_update = true;
                        // Check for special case of pvs adjacent/overlapping nbr (This is just heuristics; to ignore some of the possible wrong cases)
                        if  (pvs_supp_kpos > -1 && (spos[c_ind]<=k+pvs_supp_kpos)) { // spcl case: pvs kmer exists in the range and is overlapping/adjacent
                            if ( (r_bind-pvs_supp_r_bind) != (spos[c_ind]-pvs_supp_kpos)) { // insertion in reads between adjacent; do not support
                                should_update = false;
                            }
                        }   
                        if (should_update) {
                            pvs_supp_kpos = spos[c_ind];
                            pvs_supp_r_bind = r_bind;
                            contig.increment_support(first+c_ind);
                        }                         
                    }
                }
            }
        }
    }
}

void Alignment::update_minimisers_support (Contig& contig) {
    UINT MINIMIZER_K = Minimizer_settings.k;
    UINT MINIMIZER_W = Minimizer_settings.w;
    /* Find starting and ending region */
    // Rank(i) returns number of set nits in [0, i). Therefore look for _rb+1 so that index of starting region can be deduced by subtracting 1.
    auto first = contig._RMreg_pos(_rb+1)-1;
    auto last = contig._RMreg_pos(_re);

    auto first_windex = ((contig._is_win_even && first%2==0) || (!contig._is_win_even && first%2==1)) ? (first) : (first+1);
    auto last_windex = ((contig._is_win_even && last%2==0) || (!contig._is_win_even && last%2==1)) ? (last) : (last-1);

    // Update minimser support
    if (last_windex >= first_windex) {
        // Find minimsers
        UINT32 last_found_position = _apseq.get_seq_size() + 1; //a unique identifier for 'first minimizer'
        
        UINT32 shift = 2 * (MINIMIZER_K - 1);
        UINT32 mask = (1ULL<<2*MINIMIZER_K) - 1;
        UINT32 kmer[2] = {0,0};
        MinimizerDeque<UINT32> minimizer_window(MINIMIZER_W + 1);
        UINT32 count_not_N = 0;
        UINT32 processed_kmer = 0;
        UINT32 current_start_position = 0;
        
        //we have to make a vector of found minimizers and a map to validity, to handle duplicates
        //is there a better solution?
        std::unordered_multimap<UINT32,UINT32> found_minimizers; // maps read minimser to its position in the read
        
        for(size_t i = 0; i < _apseq.get_seq_size(); ++i) {
            BYTE c = _apseq.enc_base_at(i);
            if(c < 4) {
                ++count_not_N;                
                kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                //kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift; // reverse k-mer
                //int z = kmer[0] < kmer[1] ? 0 : 1;
                int z = 0;
                if(count_not_N >= MINIMIZER_K) {
                    while(!minimizer_window.empty() && std::get<0>(minimizer_window.back()) > kmer[z]) minimizer_window.pop_back();
                    minimizer_window.push_back(std::make_tuple(kmer[z], i));
                    while(std::get<1>(minimizer_window.front()) + MINIMIZER_W <= i) minimizer_window.pop_front();
                    ++processed_kmer;
                    if(processed_kmer >= MINIMIZER_W) {
                        current_start_position = std::get<1>(minimizer_window.front()) - MINIMIZER_K + 1;
                        if(current_start_position != last_found_position) { //first minimizer
                            found_minimizers.insert({std::get<0>(minimizer_window.front()),current_start_position});
                        }
                        last_found_position = current_start_position;
                    }                    
                }
            } else {
                count_not_N = 0;
            }
        }

        // Check minimser in each window; incrementing support if it falls within reasonable range (2*mk on either side)
        UINT16 num_cbases = _re-_rb;
        for (UINT32 i=first_windex; i<=last_windex; i=i+2) {
            UINT32 minfoidx = (contig._is_win_even) ? (i/2) : ((i-1)/2);
            auto num_minimsers = contig._minimserinfo[minfoidx]->rel_pos.size();
            // NOTE SELECT query's index is 1-based
            auto minimiser_pos = contig._SMreg_pos(i+1);         
            for (UINT32 mi=0; mi< num_minimsers; ++mi) {
                minimiser_pos+=contig._minimserinfo[minfoidx]->rel_pos[mi];
                UINT32 c_dist = minimiser_pos - _rb;
                UINT32 range_left = ((c_dist>(2*MINIMIZER_K))?(c_dist-(2*MINIMIZER_K)):(0));
                UINT32 range_right = std::min(num_cbases,(UINT16)(c_dist+(3*MINIMIZER_K)));
                if (minimiser_pos>=_rb && minimiser_pos<_re) {
                    // coverage
                    contig.increment_minimser_coverage(minfoidx,mi);
                    // support
                    auto mit = found_minimizers.equal_range(contig._minimserinfo[minfoidx]->minimisers[mi]);
                    for  (auto itr=mit.first; itr!=mit.second;++itr) {
                       
                        // WARNING: The same read-kmer can provide support to many in the contig
                        // (if it is not unique and repeat close to each other). But such cases  are rare. So ignore.
                        if (itr->second >= range_left && itr->second <= range_right) { // supports
                            contig.increment_minimser_support(minfoidx,mi);
                        }
                    }
                }
                if (minimiser_pos>=_re) {
                    break;
                }
            }
        }
    }    
}
// For adding to short windows
void Alignment::find_short_arms(const UINT k, Contig& contig) {
    /* Find the range */
    auto b_ind = contig._Rreg_pos(_rb);
    if (contig._reg_pos[_rb]==0) { // starts before
        --b_ind;
    }
    // e_ind will point to 1 past the last region falling in the read (i.e. starting pos of the next reg)
    auto e_ind = contig._Rreg_pos(_re); 
    // Ignore if the whole read falls into a single SR or window
    if (e_ind-b_ind > 1) { 
        // Find breaking-points of the read; at least 1 will exist
        std::vector<UINT32> bp = find_bp(contig._Sreg_pos,contig._reg_type,b_ind,e_ind);
        _arms.reserve(bp.size());
        // first arm may be suffix or internal or in an SR
        ArmType armtype = (contig._reg_pos[_rb]==0) ? ArmType::SUFFIX: ArmType::INTERNAL;
        if (contig._reg_type[b_ind] != RegionType::SR && contig._reg_type[b_ind] != RegionType::MSR) { // a window
            prepare_short_arm(k, b_ind, _qab, bp[0], armtype, contig);
        }
        // All others are internal arms (to be added) if in a window; Check for empty arm; If empty just incease the pointer
        UINT bp_ind = 0;
        for (auto ind=b_ind+1; ind < e_ind-1; ++ind,++bp_ind) {
            if (contig._reg_type[ind] != RegionType::SR && contig._reg_type[ind] != RegionType::MSR) { // a window
                if (bp[bp_ind+1] == bp[bp_ind]) {// Empty seq in this window
                    _arms.emplace_back(Arm(ind));
                }
                else {
                    prepare_short_arm(k, ind, bp[bp_ind], bp[bp_ind+1], ArmType::INTERNAL, contig);
                }
            }
        }
        // last window may be suffix or internal or in an SR
        armtype = (contig._reg_pos[_re]==0) ? ArmType::PREFIX: ArmType::INTERNAL;
        if (contig._reg_type[e_ind-1] != RegionType::SR && contig._reg_type[e_ind-1] != RegionType::MSR) { // a window
            prepare_short_arm(k, e_ind-1,  bp[bp_ind], _qae, armtype, contig);
        }
    }
    _arms.shrink_to_fit(); 
}

// For adding to long windows
void Alignment::find_long_arms(Contig& contig) {
    /* Find the range */
    auto b_ind = contig._pseudo_Rreg_pos(_rb);
    if (contig._pseudo_reg_pos[_rb]==0) { // starts before
        --b_ind;
    }
    // e_ind will point to 1 past the last region falling in the read (i.e. starting pos of the next reg)
    auto e_ind = contig._pseudo_Rreg_pos(_re); 
    // Ignore if the whole read falls into a single SR or window
    if (e_ind-b_ind > 1) { 
        // Find breaking-points of the read; at least 1 will exist
        std::vector<UINT32> bp = find_bp(contig._pseudo_Sreg_pos,contig._pseudo_reg_type,b_ind,e_ind);
        _arms.reserve(bp.size());
        // first arm may be suffix or internal or in an SR
        ArmType armtype = (contig._pseudo_reg_pos[_rb]==0) ? ArmType::SUFFIX: ArmType::INTERNAL;
        if (contig._pseudo_reg_type[b_ind] != RegionType::SR) { // a window (pseudo doesn't have MSR)
            _arms.emplace_back(Arm(contig._true_reg_id[b_ind], _apseq, _qab, bp[0], armtype));
        }
        // All others are internal arms (to be added) if in a window; Check for empty arm; If empty just incease the pointer
        UINT bp_ind = 0;
        for (auto ind=b_ind+1; ind < e_ind-1; ++ind,++bp_ind) {
            if (contig._pseudo_reg_type[ind] != RegionType::SR) { // a window
                if (bp[bp_ind+1] == bp[bp_ind]) {// Empty seq in this window
                    _arms.emplace_back(Arm(contig._true_reg_id[ind]));
                }
                else {
                    _arms.emplace_back(Arm(contig._true_reg_id[ind], _apseq,bp[bp_ind], bp[bp_ind+1], ArmType::INTERNAL));
                }
            }
        }
        // last window may be suffix or internal or in an SR
        armtype = (contig._pseudo_reg_pos[_re]==0) ? ArmType::PREFIX: ArmType::INTERNAL;
        if (contig._pseudo_reg_type[e_ind-1] != RegionType::SR) { // a window
            _arms.emplace_back(Arm(contig._true_reg_id[e_ind-1], _apseq, bp[bp_ind], _qae, armtype));
        }
    } 
    _arms.shrink_to_fit(); 
}

void Alignment::add_arms(const Contig& contig) {
    //std::cout << _qname <<" " << _arms.size()<<std::endl;
    for (const auto& a: _arms) {
        if (a.armtype==ArmType::PREFIX) {
            contig._pwindows[a.windex]->add_prefix(a.arm);
        }
        else if (a.armtype==ArmType::SUFFIX) {
            contig._pwindows[a.windex]->add_suffix(a.arm);
        }
        else if (a.armtype==ArmType::INTERNAL) {
            contig._pwindows[a.windex]->add_internal(a.arm);
        }
        else {
            contig._pwindows[a.windex]->add_empty();
        }
    }
    _arms.clear();
}

// beg_ind points to the starting pos of the window the read starts in; end_ind points to the starting pos of the window following the one in which read ends
std::vector<UINT32> Alignment::find_bp(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind) {    
    std::vector<UINT32> results;
    UINT32 num_bp = end_ind-beg_ind-1; // ending and beginning themselves do not fall within the read; so no bp
    results.reserve(num_bp);
    UINT32 current_reference_pos = _rb;
    UINT32 current_processed_index = beg_ind + 1;
    UINT32 next_ref_pos = Sreg_pos(current_processed_index + 1);
    UINT32 current_query_pos = 0;

    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    UINT32 len_diff;
    INT8 is_corner = 0;
    for(UINT32 i = 0; i < _num_cigar_op; i++) {        
        current_cigar = _cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        
        if(get_op == BAM_CSOFT_CLIP || get_op == BAM_CHARD_CLIP) { //softclips are already discarded
            continue;
        }
        
        if((bam_cigar_type(get_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
            }
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                current_query_pos += len_diff;
                get_oplen -= len_diff;
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                }
                else is_corner = 1;
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
                current_query_pos += get_oplen;
            }
        } else if(bam_cigar_type(get_op) & 2) { //only bit 2 set, consumes reference
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
            }
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                get_oplen -= len_diff;
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                }
                else is_corner = 1;
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
            }
        } else if(bam_cigar_type(get_op) & 1) { //only bit 1 set, consumes query
            if(is_corner) { //this is a corner, now add the position
                if(reg_type[current_processed_index-1] == RegionType::SR || reg_type[current_processed_index-1] == RegionType::MSR) { // if this is SR, we want the current_query_pos to be included in the right window
                    results.push_back(current_query_pos);
                } else { // otherwise we include it here
                    results.push_back(current_query_pos + get_oplen);
                }
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
                is_corner = 0;
            }
            current_query_pos += get_oplen;
        }
        
        if(current_processed_index == end_ind) break; //this should not be needed if end_ind is correct
    }

    return results;
}

void Alignment::prepare_short_arm(const UINT k, const UINT32 windex, const UINT32 qb, const UINT32 qe, const ArmType armtype, Contig& contig) {
    //std:: cout << "wind: qb: qe : " << windex << " " << qb << " " << qe <<" "<<std::endl;
    
    // Ignore if the arm is too short as compared to window
    auto mk = Minimizer_settings.k;
    auto curr_pos = contig._Sreg_pos(windex+1);
    auto next_pos = contig._Sreg_pos(windex+2);
    
    if ((next_pos-curr_pos) > (Arms_settings.short_arm_coef*(qe-qb))) {return;}
    RegionType wtype = contig._reg_type[windex];
    bool valid = true;
    auto q_beg = qb;
    auto q_end = qe;
    
    // Find preceding anchor SR kmer; range =/- k on either side of expected
    if ((wtype==RegionType::SWS || wtype==RegionType::SW || wtype==RegionType::SWM) && armtype!=ArmType::SUFFIX) { // preceded by SR and int or pre type
        // pa->_qab is always 0 (as clipped bases are discarded)
        if (q_beg < k) { // preceeding kmer does not exist
            valid = false;
        }
        else {
            // Get preceeding kmer ID
            UINT32 prec_SR_rank = contig._reg_info[windex-1];
            UINT64 anchor_kmer = contig._anchor_kmers[(prec_SR_rank<<1)]; // index of last kmer of SR is 2i
            // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its rightmost position
            if (!_apseq.check_kmer(anchor_kmer,k,q_beg-k)) {
                UINT32 search_start = (q_beg < (2*k)) ? (0) : (q_beg - (2*k));
                UINT32 search_end = (q_end < (q_beg+k)) ? (q_end) : (q_beg+k);
                size_t kmer_ind=0;
                if (_apseq.find_kmer(anchor_kmer,k, search_start,search_end,false,kmer_ind)) { // found
                    q_beg = kmer_ind + k;
                }
                else {valid = false;}
            }  
        }
    }
    // Find succeeding anchor SR kmer; range =/- k on either side of expected
    if ((wtype==RegionType::SWS || wtype==RegionType::WS || wtype==RegionType::MWS) && armtype!=ArmType::PREFIX) { // succeeded by SR and int or suff type
        // pa->_qab is always 1 (as clipped bases rae discarded)
        if (q_end + k > _qae) { // succeeding kmer does not exist
            valid = false;
        }
        else {
            // Get succeeding kmer ID
            UINT32 succ_SR_rank = contig._reg_info[windex+1];
            UINT64 anchor_kmer = contig._anchor_kmers[(succ_SR_rank<<1)-1]; // index of first kmer of SR is 2i-1
            // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its leftmost position
            if (!_apseq.check_kmer(anchor_kmer,k,q_end)) {
                UINT32 search_start = (q_end < (q_beg+k)) ? (q_beg) : (q_end-k);
                UINT32 search_end = std::min(_qae,q_end+(2*k));
                size_t kmer_ind=0;
                if (_apseq.find_kmer(anchor_kmer, k,search_start,search_end,true,kmer_ind)) { // found
                    q_end = kmer_ind;
                }
                else {valid = false;}
            }
        }
    }
    // Find preceding anchor minimiser kmer; range =/- 2*mk on either side of expected
    if ((wtype==RegionType::MWM || wtype==RegionType::MW || wtype==RegionType::MWS) && armtype!=ArmType::SUFFIX) { // preceded by MINI and int or pre type
        // pa->_qab is always 0 (as clipped bases are discarded)
        if (q_beg < mk) { // preceeding kmer does not exist
            valid = false;
        }
        else {
            // Get preceeding kmer ID
            UINT32 prec_MIN = contig._reg_info[windex-1];
            // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its rightmost position
            if (!_apseq.check_kmer(prec_MIN,mk,q_beg-mk)) {
                UINT32 search_start = (q_beg < (3*mk)) ? (0) : (q_beg - (3*mk));
                UINT32 search_end = (q_end < (q_beg+(2*mk))) ? (q_end) : (q_beg+(2*mk));
                size_t kmer_ind=0;
                if (_apseq.find_kmer(prec_MIN,mk, search_start,search_end,false,kmer_ind)) { // found
                    q_beg = kmer_ind + mk;
                }
                else {valid = false;}
            }
        }
    }
    // Find succeeding minimiser kmer; range =/- 2*mk on either side of expected
    if ((wtype==RegionType::MWM || wtype==RegionType::WM || wtype==RegionType::SWM) && armtype!=ArmType::PREFIX) { // succeeded by MINI and int or suff type
        // pa->_qab is always 1 (as clipped bases rae discarded)
        if (q_end + mk > _qae) { // succeeding kmer does not exist
            valid = false;
        }
        else {
            // Get succeeding kmer ID
            UINT32 succ_MIN = contig._reg_info[windex+1];
            // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its leftmost position
            if (!_apseq.check_kmer(succ_MIN,mk,q_end)) {
                UINT32 search_start = (q_end < (q_beg+(2*mk))) ? (q_beg) : (q_end-(2*mk));
                UINT32 search_end = std::min(_qae,q_end+(3*mk));
                size_t kmer_ind=0;
                if (_apseq.find_kmer(succ_MIN, mk,search_start,search_end,true,kmer_ind)) { // found
                    q_end = kmer_ind;
                }
                else {valid = false;}
            }
        }
    }
    if (valid && q_beg < q_end) {
        _arms.emplace_back(Arm(windex, _apseq, q_beg,q_end,armtype));
    }
}


void Alignment::initialise_pos(const bam1_t *hts_align) {
    //_qname = std::string(bam_get_qname(hts_align), hts_align->core.l_qname);
    _rb = hts_align->core.pos;
    _qab = 0;
    // Get ending pos and update _qab (if soft-clipped)
    
    UINT32 curr_qp = _qab;
    UINT32 curr_rp = _rb;
    bool clip_before = true;
    UINT32 * pcigar = bam_get_cigar(hts_align);
    _num_cigar_op = hts_align->core.n_cigar;
    UINT32 count_clip_end = 0;
    for(UINT32 i = 0; i < _num_cigar_op; ++i) {
        UINT32 current_cigar = pcigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);        
        if(clip_before) {
            if(cigar_op == BAM_CSOFT_CLIP) {
                _qab += cigar_oplen;
            } else if(cigar_op != BAM_CHARD_CLIP) {
                clip_before = false;
            }
        }        
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            curr_rp += cigar_oplen;
            curr_qp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            curr_rp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            if(!clip_before && cigar_op == BAM_CSOFT_CLIP) {count_clip_end += cigar_oplen;}
            curr_qp += cigar_oplen;
        }
    }
    _qae = curr_qp - count_clip_end;
    _re = curr_rp; 
}

void Alignment::copy_data(const bam1_t *hts_align) {
    UINT32 * pcigar = bam_get_cigar(hts_align);
    // Get seq
    UINT32 qlen = _qae-_qab;
    UINT32 offset = _qab;
    PackedSeq<2> ps (qlen,offset,bam_get_seq(hts_align));
    if (ps.is_valid()) {
        _apseq = std::move(ps);
        // Reset query indices in accordance with only the stored (aligned) part; Clipped part already discarded
        _qae = _qae-_qab;
        _qab = 0;
        
        // Save cigar (without changing HTS format)
        // Cigar data is encoded 4 bytes per CIGAR operation.
        UINT64 cigar_data_size = _num_cigar_op << 2; // multiply by 4
        _cigar.reserve(cigar_data_size);
        _cigar.insert(_cigar.begin(),pcigar,pcigar+cigar_data_size); // copy not move; else htslib will destroy it  
    }
    else {is_valid=false;}
      
}


} // namespace hypo
