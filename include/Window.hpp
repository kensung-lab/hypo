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
/** Class Window.
 * It contains the functionality for POA.
 */
#pragma once
#include "spoa/spoa.hpp"
#include "globalDefs.hpp"
#include "PackedSeq.hpp"
#include "Filter.hpp"

namespace hypo
{
const std::string cMarkerLetters = "JO";
const std::string cHead = "J";
const std::string cTail = "O";
const UINT cMarkerLen = cHead.size();

enum class WindowType: UINT8 {
    SHORT,
    LONG
};


class Window {
public:
    Window(): _wtype(WindowType::SHORT), 
    _num_internal(0), _num_pre(0), _num_suf(0), _num_empty(0),
    _longest_pre_len(0),_longest_suf_len(0), _filter() {}

    Window(const PackedSeq<4> & ps, const size_t left_ind, const size_t right_ind, WindowType wt): _wtype(wt), 
    _num_internal(0), _num_pre(0), _num_suf(0), _num_empty(0), _longest_pre_len(0),_longest_suf_len(0), _filter(), _draft(ps,left_ind,right_ind)
    {
        if (wt==WindowType::LONG) {
            _filter.initialise(ps.unpack(left_ind,right_ind));
        }
    }

    Window(const Window &) = delete;
    Window &operator=(const Window &) = delete;
	Window(Window&&) = delete;
    Window& operator=(Window&&) = delete;
    ~Window() = default;

    static void prepare_for_poa(const ScoreParams& sp, const UINT32 num_threads);
    void generate_consensus(const UINT32 engine_idx);
    inline std::string get_consensus() const {return _consensus;}
    inline size_t get_window_len() const {return _draft.get_seq_size();}

    inline void add_prefix(const PackedSeq<2>& ps) { 
        bool should_add = true;
        if (_wtype==WindowType::LONG) {
            should_add = _filter.is_good(ps.unpack());
        }
        if (should_add) {
            UINT arm_len = ps.get_seq_size();         
            ++_num_pre;
            if (arm_len > _longest_pre_len) {_longest_pre_len= arm_len;}
            _pre_arms.emplace_back(std::move(ps));
        }       
    }

    inline void add_suffix(const PackedSeq<2>& ps) { 
        bool should_add = true;
        if (_wtype==WindowType::LONG) {
            should_add = _filter.is_good(ps.unpack());
        }
        if (should_add) {
            UINT arm_len = ps.get_seq_size();
            ++_num_suf;
            if (arm_len > _longest_suf_len) {_longest_suf_len= arm_len;}
            _suf_arms.emplace_back(std::move(ps)); 
        }             
    }

    inline void add_internal(const PackedSeq<2>& ps) { 
        bool should_add = true;
        if (_wtype==WindowType::LONG) {
            should_add = _filter.is_good(ps.unpack());
        }
        if (should_add) {
            ++_num_internal;
            _internal_arms.emplace_back(std::move(ps)); 
        }        
    }

    inline void add_empty() {++_num_empty;}

    inline UINT32 get_num_pre() const {return _num_pre;}
    inline UINT32 get_num_suf() const {return _num_suf;}
    inline UINT32 get_num_internal() const {return _num_internal+_num_empty;}
    inline UINT32 get_num_total() const {return _num_internal+_num_empty+_num_pre+_num_suf;}
    inline UINT32 get_maxlen_pre() const {return _longest_pre_len;}
    inline UINT32 get_maxlen_suf() const {return _longest_suf_len;}
    inline void clear_pre_suf() {
        _num_pre = 0;
        _num_suf = 0;
        _pre_arms.clear();
        _suf_arms.clear();
        _pre_arms.shrink_to_fit();
        _suf_arms.shrink_to_fit();
    }

    friend std::ostream &operator<<(std::ostream &, const Window &);

private:
    WindowType _wtype;
    UINT32 _num_internal;
    UINT32 _num_pre;
    UINT32 _num_suf;
    UINT32 _num_empty;
    UINT32 _longest_pre_len;
    UINT32 _longest_suf_len;
    PackedSeq<4> _draft;
    std::vector<PackedSeq<2>> _internal_arms;
    std::vector<PackedSeq<2>> _pre_arms;
    std::vector<PackedSeq<2>> _suf_arms;
    std::string _consensus;
    Filter _filter; // only for long windows

    static std::vector<std::shared_ptr<spoa::AlignmentEngine>> _alignment_engines;
    static std::vector<std::shared_ptr<spoa::AlignmentEngine>> _alignment_engines_long;
    static const float _cThresh;

    void generate_consensus_short(const UINT32 engine_idx);
    void generate_consensus_long(const UINT32 engine_idx, bool initial);
    inline void set_consensus(std::string const &con){_consensus = std::move(con);}
    inline void set_marked_consensus(std::string const &con){_consensus.assign(con.begin()+cMarkerLen,con.end()-cMarkerLen);}
    std::string curate(const std::string &con, const std::vector<UINT32> &dst);
}; // Window
} // namespace hypo