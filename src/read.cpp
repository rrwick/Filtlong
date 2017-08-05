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
#include <math.h>
#include <limits>
#include <sstream>
#include <iomanip>

#include "read.h"

Read::Read(char * name, char * seq, char * qscores, int length, Kmers * kmers, int window_size) {
    m_name = name;
    m_length = length;

    std::vector<double> qualities;

    // If reference k-mers aren't available, use the qscores to get the qualities.
    if (kmers->empty()) {
        qualities.reserve(length);
        for (int i = 0; i < length; ++i)
            qualities.push_back(qscore_to_quality(qscores[i]));
    }

    // If there are reference k-mers, use them for the qualities. A base is considered to have a quality of 1 if it
    // is in any present 16-mer, 0 if it is not.
    else {
        qualities.resize(length, 0.0);
        uint32_t kmer;
        if (length >= 16) {
            for (int i = 15; i < length; ++i) {
                if (i == 15)
                    kmer = kmers->starting_kmer_to_bits_forward(seq);
                else {
                    kmer <<= 2;
                    kmer |= kmers->base_to_bits_forward(seq[i]);
                }
                if (kmers->is_kmer_present(kmer)) {
                    for (int j = i - 15; j <= i; ++j)
                        qualities[j] = 1.0;
                }
            }
        }
    }

    m_mean_quality = get_mean_quality(qualities);
    m_window_quality = get_window_quality(qualities, window_size);

    m_length_score = 0.0;
    m_mean_quality_score = 0.0;
    m_window_quality_score = 0.0;
    m_final_score = 0.0;
    m_passed = true;
}


std::string double_to_string(double a) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << a;
    return ss.str();
}


void Read::print_table_row() {
    std::cerr << m_name << "\t";
    std::cerr << m_length << "\t";
    std::cerr << double_to_string(m_mean_quality) << "\t";
    std::cerr << double_to_string(m_window_quality) << "\t";
    std::cerr << double_to_string(m_length_score) << "\t";
    std::cerr << double_to_string(m_mean_quality_score) << "\t";
    std::cerr << double_to_string(m_window_quality_score) << "\t";
    std::cerr << double_to_string(m_final_score) << "\t";
    std::cerr << m_passed << "\t";
    std::cerr << "\n";
}


double Read::get_mean_quality(std::vector<double> & qualities) {
    double sum = 0.0;
    for (auto q : qualities)
        sum += q;
    return 100.0 * sum / qualities.size();
}


double Read::get_window_quality(std::vector<double> & qualities, int window_size) {
    if (qualities.size() <= window_size)
        return get_mean_quality(qualities);

    double sum = 0.0;
    for (int i = 0; i < window_size; ++i)
        sum += qualities[i];
    double window_quality = sum / window_size;
    double min_window_quality = window_quality;

    for (int j = window_size; j < qualities.size(); ++j) {
        int i = j - window_size;
        window_quality -= qualities[i] / window_size;
        window_quality += qualities[j] / window_size;
        if (window_quality < min_window_quality)
            min_window_quality = window_quality;
    }
    if (min_window_quality < 0.5 / window_size)
        min_window_quality = 0.0;
    return 100.0 * min_window_quality;
}


double Read::qscore_to_quality(char qscore) {
    int q = qscore - 33;
    return 1.0 - pow(10.0, -q / 10.0);
}


void print_read_table_header() {
    std::cerr << "Read name" << "\t";
    std::cerr << "Length" << "\t";
    std::cerr << "Mean quality" << "\t";
    std::cerr << "Window quality" << "\t";
    std::cerr << "Length score" << "\t";
    std::cerr << "Mean quality score" << "\t";
    std::cerr << "Window quality score" << "\t";
    std::cerr << "Final score" << "\t";
    std::cerr << "Passed";
    std::cerr << "\n";
}