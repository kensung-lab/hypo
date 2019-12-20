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
/** Class MinimzerDeque.
 * Double ended queue used for minimizer finding.
 */
#pragma once

#include <unordered_map>
#include "globalDefs.hpp"

namespace hypo
{
template <typename T>
class MinimizerDeque {
           
    public:
        MinimizerDeque(UINT c) {
            if (c > 0) {
                _capacity = c;
                _front_idx = 0;
                _back_idx = c - 1;
                _count = 0;
                _content.resize(c);
            }            
        }
        
        inline std::tuple<T, UINT32> & front() {
            return _content[_front_idx];
        }
        
        inline std::tuple<T, UINT32> & back() {
            return _content[_back_idx];
        }
        
        inline UINT size() {
            return _count;
        }
        
        inline void push_back(std::tuple<T, UINT32> x) {
            ++_count;
            _back_idx = (_back_idx + 1) % _capacity;
            _content[_back_idx] = x;
        }
        
        inline void pop_back() {
            --_count;
            if(_back_idx == 0) {_back_idx = _capacity-1;}
            else {--_back_idx;}
        }
        
        inline void push_front(std::tuple<T, UINT32> x) {
            ++_count;
            if(_front_idx == 0) {_front_idx = _capacity-1;}
            else {--_front_idx;}
            _content[_front_idx] = x;
        }
        
        inline void pop_front() {
            _front_idx = (_front_idx + 1) % _capacity;
            --_count;
        }
        
        inline bool empty() {
            return _count == 0;
        }
        
        inline void clear() {
            _front_idx = 0; _back_idx = _capacity-1; _count = 0;
        }

    private:
        std::vector<std::tuple<T, UINT32> > _content;
        UINT _front_idx;
        UINT _back_idx;
        UINT _capacity;
        UINT _count;
};


}