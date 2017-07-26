//Copyright 2017 Ryan Wick

//This file is part of LongQC

//LongQC is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//LongQC is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with LongQC.  If not, see <http://www.gnu.org/licenses/>.


#include <iostream>
#include "read.h"

Read::Read(char * name, char * seq, char * quals, int length, Kmers * kmers, int window_size) {
    m_name = name;
    m_length = length;

    // If reference k-mers aren't available, use the read quality scores.
    if (kmers) {
        m_mean_quality = get_mean_quality(quals);
        m_window_quality = get_window_quality(quals, window_size);
    }

    // If there are reference k-mers, use them instead of read quality scores.
    else {

    }
}


void Read::print_table_row() {
    std::cerr << m_name << "\t";
    std::cerr << m_length << "\t";
    // TO DO: rest of table
    std::cerr << "\n";
}


double Read::get_mean_quality(char * quals) {
    return 0.0; // TEMP
}


double Read::get_window_quality(char * quals, int window_size) {
    return 0.0; // TEMP
}


void print_read_table_header() {
    std::cerr << "Read name" << "\t";
    std::cerr << "Length" << "\t";
    std::cerr << "Mean quality" << "\t";
    std::cerr << "Window quality" << "\t";
    std::cerr << "Length score" << "\t";
    std::cerr << "Mean quality score" << "\t";
    std::cerr << "Window quality score" << "\t";
    std::cerr << "Total quality" << "\t";
    std::cerr << "Passed";
    std::cerr << "\n";
}