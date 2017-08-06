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
#include <string>
#include <iomanip>

#include "read.h"

Read::Read(char * name, char * seq, char * qscores, int length, Kmers * kmers, Arguments * args) {
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
    m_window_quality = get_window_quality(qualities, args->window_size);

    m_length_score = get_length_score();
    m_final_score = get_final_score(args->length_weight, args->mean_q_weight, args->window_q_weight);

    // See if the read failed any of the hard cut-offs.
    m_passed = true;
    if (m_length < args->min_length)
        m_passed = false;
    else if (m_mean_quality < args->min_mean_q)
        m_passed = false;
    else if (m_window_quality < args->min_window_q)
        m_passed = false;

    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
    // TO DO: trimming/splitting
}


std::string double_to_string(double a) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << a;
    std::string s = ss.str();
    if (s.size() < 5)
        return std::string(5 - s.size(), ' ') + ss.str();
    else
        return ss.str();
}


std::string pad(int num, const size_t width)
{
    std::string s = std::to_string(num);
    if (width > s.size())
        return s + std::string(width - s.size(), ' ');
    else
        return s;
}


void Read::print_verbose_read_info() {
    std::cerr << "\n" << m_name << "\n";

    std::cerr << "            length = " << pad(m_length, 11);
    std::cerr << "length score = " << double_to_string(m_length_score) << "\n";

    std::cerr << "      mean quality = " << double_to_string(m_mean_quality);
    std::cerr << "    window quality = " << double_to_string(m_window_quality);

    std::cerr << "   final score = " << double_to_string(m_final_score) << "\n";
}


double Read::get_mean_quality(std::vector<double> & qualities) {
    double sum = 0.0;
    for (auto q : qualities)
        sum += q;
    return 100.0 * sum / qualities.size();
}


double Read::get_window_quality(std::vector<double> & qualities, size_t window_size) {
    if (qualities.size() <= window_size)
        return get_mean_quality(qualities);

    double sum = 0.0;
    for (size_t i = 0; i < window_size; ++i)
        sum += qualities[i];
    double window_quality = sum / window_size;
    double min_window_quality = window_quality;

    for (size_t j = window_size; j < qualities.size(); ++j) {
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

// At the moment, the half-score length is hard-coded to 5 kbp. Maybe this should be adjustable via a setting?
// https://www.desmos.com/calculator
// y=100\left(1+\frac{-a}{x+a}\right)
double Read::get_length_score() {
    double half_length_score = 5000.0;
    return 100.0 * (1.0 + (-half_length_score / (m_length + half_length_score)));
}


// The final score is a weighted geometric mean of the length score and the mean quality. It is then scaled down using
// the window quality.
double Read::get_final_score(double length_weight, double mean_q_weight, double window_q_weight) {

    // First get the weighted geometric mean of the length score and the mean quality.
    double product = pow(m_length_score, length_weight) * pow(m_mean_quality, mean_q_weight);
    double total_weight = length_weight + mean_q_weight;
    double final_score = pow(product, 1.0 / total_weight);

    // Now scale that down using the ratio of window quality to mean quality.
    double scaling_factor = m_window_quality / m_mean_quality;
    if (scaling_factor > 1.0)
        scaling_factor = 1.0;
    total_weight = length_weight + mean_q_weight + window_q_weight;
    double window_weight_fraction = window_q_weight / total_weight;
    double non_window_weight_fraction = 1.0 - window_weight_fraction;
    scaling_factor = (1.0 * non_window_weight_fraction) + (scaling_factor * window_weight_fraction);
    return final_score * scaling_factor;
}


double Read::qscore_to_quality(char qscore) {
    int q = qscore - 33;
    return 1.0 - pow(10.0, -q / 10.0);
}
