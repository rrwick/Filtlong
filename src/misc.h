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

#ifndef MISC_H
#define MISC_H

#include <string>

std::string double_to_string(double n);
std::string int_to_string(long long n);
void print_hash_progress(std::string filename, long long base_count);
void print_read_score_progress(int read_count, long long base_count);


#endif // MISC_H
