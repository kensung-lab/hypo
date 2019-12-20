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
/** Class Filter.
 * It filters out long arms that are too different from the draft in a window.
 */
#pragma once

#include <unordered_map>
#include "globalDefs.hpp"
#include "MinimizerDeque.hpp"
namespace hypo
{

class Filter {
            
    public:
        void initialise(const std::string & draft) {
            UINT64 shift = 2 * (_cMinimizer_k - 1);
            UINT64 mask = (1ULL<<2*_cMinimizer_k) - 1;
            UINT64 kmer[2] = {0,0};
            MinimizerDeque<UINT64> minimizer_window(_cMinimizer_w + 1);
            UINT count_not_N = 0;
            UINT processed_kmer = 0;
            for(size_t i = 0; i < draft.size(); ++i) {
                BYTE c = cNt4Table[(BYTE)draft[i]];
                if(c < 4) {
                    ++count_not_N;
                    
                    kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                    kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift; // reverse k-mer
                    int z = kmer[0] < kmer[1] ? 0 : 1;
                    if(count_not_N >= _cMinimizer_k) {
                        while(!minimizer_window.empty() && std::get<0>(minimizer_window.back()) > kmer[z]) {minimizer_window.pop_back();}
                        minimizer_window.push_back(std::make_tuple(kmer[z], i));
                        while(std::get<1>(minimizer_window.front()) + _cMinimizer_w <= i) {minimizer_window.pop_front();}
                        ++processed_kmer;
                        if(processed_kmer >= _cMinimizer_w) {
                            _draft_minimizers.insert(std::get<0>(minimizer_window.front()));
                        }
                    }
                } 
                else {
                    count_not_N = 0;
                }
            }
        }
        
        bool is_good(const std::string & arm) {
            UINT64 shift = 2 * (_cMinimizer_k - 1);
            UINT64 mask = (1ULL<<2*_cMinimizer_k) - 1;
            UINT64 kmer[2] = {0,0};
            MinimizerDeque<UINT64> minimizer_window(_cMinimizer_w + 1);
            UINT count_not_N = 0;
            UINT processed_kmer = 0;
            std::vector<std::tuple<UINT64, UINT32> > found_minimizers;
            for(size_t i = 0; i < arm.size(); ++i) {
                BYTE c = cNt4Table[(BYTE)arm[i]];
                if(c < 4) {
                    ++count_not_N;
                    
                    kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                    kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift; // reverse k-mer
                    int z = kmer[0] < kmer[1] ? 0 : 1;
                    if(count_not_N >= _cMinimizer_k) {
                        while(!minimizer_window.empty() && std::get<0>(minimizer_window.back()) > kmer[z]) {minimizer_window.pop_back();}
                        minimizer_window.push_back(std::make_tuple(kmer[z], i));
                        while(std::get<1>(minimizer_window.front()) + _cMinimizer_w <= i) {minimizer_window.pop_front();}
                        
                        ++processed_kmer;
                        if(processed_kmer >= _cMinimizer_w) {
                            if(found_minimizers.size() == 0 || std::get<1>(found_minimizers[found_minimizers.size()-1]) != std::get<1>(minimizer_window.front())) {
                                found_minimizers.push_back(minimizer_window.front());
                            }
                        }
                    }
                } 
                else {
                    count_not_N = 0;
                }
            }
            uint32_t count_found_minimizer = 0;
            for(size_t i = 0; i < found_minimizers.size(); i++) {
                if(_draft_minimizers.find(std::get<0>(found_minimizers[i])) != _draft_minimizers.end()) {++count_found_minimizer;}
            }
            return (count_found_minimizer * _cThreshold_bp >= arm.size());
        }

    private:
        std::unordered_set<UINT64> _draft_minimizers;
        static const UINT _cMinimizer_k = 10;  //k-length minimizer
        static const UINT _cMinimizer_w = 10;  //window size of minimizer
        static const UINT _cThreshold_bp = 50; //1 minimizer every 50 bp
        
};
}