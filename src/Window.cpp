/*
 * 
 * Copyright (c) 2019, Ritu Kundu  and Joshua Casey
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
/** Defines the class Window.
 * It contains the functionality for POA.
 */
#include <cmath>
#include <unordered_map>
#include "Window.hpp"

namespace hypo
{
    const float Window::_cThresh = 0.4;
    std::vector<std::shared_ptr<spoa::AlignmentEngine>> Window::_alignment_engines;
    std::vector<std::shared_ptr<spoa::AlignmentEngine>> Window::_alignment_engines_long;
    void Window::prepare_for_poa(const ScoreParams& sp, const UINT32 num_threads) {
        for (UINT32 i = 0; i < num_threads; ++i){
            _alignment_engines.emplace_back(spoa::createAlignmentEngine(
                spoa::AlignmentType::kNW, sp.sr_match_score,
                sp.sr_misMatch_score, sp.sr_gap_penalty));
            _alignment_engines_long.emplace_back(spoa::createAlignmentEngine(
                spoa::AlignmentType::kNW, sp.lr_match_score,
                sp.lr_misMatch_score, sp.lr_gap_penalty));
        }
        _alignment_engines.shrink_to_fit();
        _alignment_engines_long.shrink_to_fit();
    } 

void Window::generate_consensus(const UINT32 engine_idx) {
    auto wt = _wtype;
    auto num_non_empty_arms = _num_internal+_num_pre+_num_suf;
    if (_num_empty > num_non_empty_arms) {
        set_consensus(""); // empty sequence is consensus
    }
    else if (num_non_empty_arms >= 2) { // 1 is draft; there should be at least 2 arms to do poa
        if (wt == WindowType::SHORT) {
            generate_consensus_short(engine_idx);
        }
        else {
            generate_consensus_long(engine_idx,true);
        }
    }
    else {
        set_consensus(_draft.unpack());
    }  
}

std::ostream &operator<<(std::ostream &os, const Window &wnd) {
        // Write numbers
        os << wnd._num_internal << "\t" << wnd._num_pre << "\t" << wnd._num_suf << "\t" << wnd._num_empty << std::endl;
        // Write draft sequence
        os << "++\t" << wnd._draft.unpack() << std::endl;

        // Write consensus sequence
        os << "++\t" << wnd._consensus << std::endl;
        // Write internal arms
        for (UINT32 i=0; i < wnd._num_internal; ++i) {
            os << wnd._internal_arms[i].unpack() << std::endl;
        }
        // Write pref arms
        for (UINT32 i=0; i < wnd._num_pre; ++i) {
            os << wnd._pre_arms[i].unpack() << std::endl;
        }
        // Write suff arms
        for (UINT32 i=0; i < wnd._num_suf; ++i) {
            os << wnd._suf_arms[i].unpack() << std::endl;
        }
        return os;
    }


void Window::generate_consensus_short(const UINT32 engine_idx) {
    std::string consensus = "";
    auto graph = spoa::createGraph();
    bool arms_added = false;
    
    // internal
    // Avoid draft if an internal arm exists
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kNW);
    if (_internal_arms.size()==0) { // use draft only when no internal arm
        // draft will never have zero size
        //std::string unpacked_str(std::move(_internal_arms[i].unpack()));
        std::string modified_str = cHead + _draft.unpack() + cTail;
        auto alignment = _alignment_engines[engine_idx]->align(modified_str, graph);
        graph->add_alignment(alignment, modified_str);
    }
    for (UINT i =0; i <  _internal_arms.size(); ++i) {
        if (_internal_arms[i].get_seq_size() > 0) {                
            //std::string unpacked_str(std::move(_internal_arms[i].unpack()));
            std::string modified_str = cHead + _internal_arms[i].unpack() + cTail;
            arms_added = true;
            auto alignment = _alignment_engines[engine_idx]->align(modified_str, graph);
            graph->add_alignment(alignment, modified_str);
        }
    }
    // prefix (add in reverse order (Since bam is sorted. last one should be the longest)
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kLOV);
    for (auto it = _pre_arms.rbegin(); it!=_pre_arms.rend(); ++it) {
        if (it->get_seq_size() > 0) {
        //std::string unpacked_str(std::move(it.unpack()));
        std::string modified_str = cHead + it->unpack();
        arms_added = true;
        auto alignment = _alignment_engines[engine_idx]->align(modified_str, graph);
        graph->add_alignment(alignment, modified_str);
        }
    }
    //suffix
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kROV);
    for (auto it = _suf_arms.begin(); it!=_suf_arms.end(); ++it) {
        if (it->get_seq_size() > 0) {
        //std::string unpacked_str(std::move(it.unpack()));
        std::string modified_str = it->unpack() + cTail;
        arms_added = true;
        auto alignment = _alignment_engines[engine_idx]->align(modified_str, graph);
        graph->add_alignment(alignment, modified_str);
        }
    }
    ///*************** FOR DEBUGGING ***************
    #ifdef DEBUG3
    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa, true);
    auto seq_num = 0;
    std::cout << "================================="<< seg.get_id() <<std::endl;
    for (const auto &it : msa){
        std::cout << seq_num<< "\t" << it.c_str() << std::endl;
        ++seq_num;
    }
    #endif
    //*/ // /********************************************/
    if (arms_added)
    { 
        consensus = graph->generate_consensus();
        set_marked_consensus(consensus);   
    }
    else {
        set_consensus(_draft.unpack());
    }
    
}

void Window::generate_consensus_long(const UINT32 engine_idx, bool initial) {
    std::string consensus = "";
    auto graph = spoa::createGraph();
    bool arms_added = false;
    auto num_non_empty_arms = _num_internal+_num_pre+_num_suf;
    
    std::vector<std::vector<UINT32>> counts;
    counts.reserve(num_non_empty_arms);
    
    // internal
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kNW);
    if (!initial) { // ignore draft and add consensus, if re-round of polishing
        if (_consensus.size() > 0) {
            auto alignment = _alignment_engines_long[engine_idx]->align(_consensus, graph);
            graph->add_alignment(alignment, _consensus);
        }

    }
    else { //  add draft
        std::string unpacked_str(std::move(_draft.unpack()));
        // draft's size will never be zero
        auto alignment = _alignment_engines_long[engine_idx]->align(unpacked_str, graph);
        graph->add_alignment(alignment, unpacked_str);
    }
    for ( UINT i=0; i <  _internal_arms.size(); ++i) {
        std::string unpacked_str(std::move(_internal_arms[i].unpack()));
        if (unpacked_str.size() > 0) {
            arms_added = true;
            auto alignment = _alignment_engines_long[engine_idx]->align(unpacked_str, graph);
            graph->add_alignment(alignment, unpacked_str);
        }
    }
    // prefix
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kLOV);
    for (const auto &it : _pre_arms) {
        std::string unpacked_str(std::move(it.unpack()));
        if (unpacked_str.size() > 0) {
        arms_added = true;
        auto alignment = _alignment_engines_long[engine_idx]->align(unpacked_str, graph);
        graph->add_alignment(alignment, unpacked_str);
        }
    }
    //suffix
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kROV);
    for (const auto &it : _suf_arms) {
        std::string unpacked_str(std::move(it.unpack()));
        if (unpacked_str.size() > 0) {
        arms_added = true;
        auto alignment = _alignment_engines_long[engine_idx]->align(unpacked_str, graph);
        graph->add_alignment(alignment, unpacked_str);
        }
    }
    ///*************** FOR DEBUGGING ***************
    #ifdef DEBUG3
    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa, true);
    auto seq_num = 0;
    std::cout << "================================="<< seg.get_id() <<std::endl;
    for (const auto &it : msa){
        std::cout << seq_num<< "\t" << it.c_str() << std::endl;
        ++seq_num;
    }
    #endif
    //*/ // /********************************************/
    if (arms_added)
    {
        std::vector<UINT32> dst;
        consensus = graph->generate_consensus_custom(dst);
        set_consensus(curate(consensus, dst));  
        //std::vector<std::string> msa;
        //graph->generate_multiple_sequence_alignment(msa, true);
        //consensus = graph->generate_consensus();
        //set_consensus(expand(consensus, dst, msa, counts),false); 
        if (initial) { // call second round of polishing
            generate_consensus_long(engine_idx,false);
        }           
    }
    else {
        set_consensus(_draft.unpack());
    } 
}


std::string Window::curate(const std::string &con, const std::vector<UINT32> &dst) {
    //auto num_arms = _num_internal+_num_pre+_num_suf+_num_empty;
    auto conssize = con.size();
    std::string curated_consensus;
    curated_consensus.reserve(conssize);

    UINT cov_thres = (UINT)std::floor(_num_internal * _cThresh);

    for (UINT i = 0; i < conssize; ++i) {
        if (dst[i] >= cov_thres) {
            curated_consensus.push_back(con[i]);
        }
    }
    curated_consensus.shrink_to_fit();
    return curated_consensus;
}
    
} // namespace hypo
