/*
 * 
 * Copyright (c) 2019, Ritu Kundu
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
#include <cmath>
#include <iostream>
#include <string>
#include "slog/Monitor.hpp"
#include "memory_usage.h"

namespace slog
{
Monitor::Monitor():
_start_tp(), _initial_tp(){}

Monitor::~Monitor() {}

void Monitor::start() {
    auto now = std::chrono::steady_clock::now();
    if (_initial_tp == std::chrono::time_point<std::chrono::steady_clock>()) { // first call
        _initial_tp = now;
    }
    _start_tp = now;
}

std::string Monitor::stop(const std::string& msg) {
    auto end = std::chrono::steady_clock::now();
    double elapsed_secs = std::chrono::duration_cast<std::chrono::duration<double> >(end - _start_tp).count();
    _start_tp = end;
    size_t currentRSS = std::ceil(getCurrentRSS()/MBYTES);
    size_t peakRSS = std::ceil(getPeakRSS()/MBYTES);
    std::cout << "RESOURCES (" << msg << "): TIME= " << elapsed_secs << " sec; PEAK RSS (so far)= " << peakRSS << "MB; CURRENT RSS (so far)= " << currentRSS << "MB."<<std::endl;
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return s;
}

std::string Monitor::total(const std::string& msg) {
    auto end = std::chrono::steady_clock::now();
    double elapsed_secs = std::chrono::duration_cast<std::chrono::duration<double> >(end - _initial_tp).count();
    _initial_tp = std::chrono::time_point<std::chrono::steady_clock>();
    size_t currentRSS = std::ceil(getCurrentRSS()/MBYTES);
    size_t peakRSS = std::ceil(getPeakRSS()/MBYTES);
    std::cout << "RESOURCES (" << msg << "): TIME= " << elapsed_secs << " sec; PEAK RSS (so far)= " << peakRSS << "MB; CURRENT RSS (so far)= " << currentRSS << "MB."<<std::endl;
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return s;
}

} // namespace slog


