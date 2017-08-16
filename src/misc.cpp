// Copyright 2017 Ryan Wick

// This file is part of Filtlong

// Filtlong is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
// version.

// Filtlong is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.

// You should have received a copy of the GNU General Public License along with Filtlong.  If not, see
// <http://www.gnu.org/licenses/>.


#include "misc.h"

#include <iostream>
#include <sstream>
#include <iomanip>


std::string double_to_string(double n) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << n;
    std::string s = ss.str();
    if (s.size() < 5)
        return std::string(5 - s.size(), ' ') + ss.str();
    else
        return ss.str();
}


std::string int_to_string(long long n) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << n;
    return ss.str();
}

void print_hash_progress(std::string filename, long long base_count) {
    std::cerr << "\r  " << filename << " (" << int_to_string(base_count) << " bp)";
}


void print_read_score_progress(int read_count, long long base_count) {
    std::cerr << "\r  " << int_to_string(read_count) << " reads (" << int_to_string(base_count) << " bp)";
}
