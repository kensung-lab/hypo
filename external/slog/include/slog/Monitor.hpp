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
/** Defines the class Logger.
 * 
 */
#pragma once
#include <chrono>

using TP = std::chrono::time_point<std::chrono::steady_clock>;

namespace slog
{
using UINT = unsigned int;
constexpr UINT MBYTES = 1024 * 1024;
class Monitor
{
public:
	
	/*!
     * @brief Sets up the monitor to print on stderr 
	 * NOT THREAD SAFE
     */
	Monitor();


	Monitor(const Monitor&) = default;
    Monitor& operator=(const Monitor&) = default;

    Monitor(Monitor&&) = default;
    Monitor& operator=(Monitor&&) = default;

    ~Monitor();
	
	void start();
	/*!
     * @brief Stops the timer and resets it.
	 * Prints the elapsed time starting from the previous call of the function measure_start.
     */
	std::string stop(const std::string& msg="");
	/*!
     * @brief Prints the elapsed time starting from the first call of the function measure_start.
     */
	std::string total(const std::string& msg="Overall");

private:
    TP _start_tp;
	TP _initial_tp;

}; // class Logger

} // namespace slog
