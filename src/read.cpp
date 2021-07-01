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


#include <iostream>
#include <math.h>
#include <limits>
#include <string>

#include "read.h"
#include "misc.h"

Read::Read(std::string name, char * seq, char * qscores, int length, Kmers * kmers, Arguments * args) {
    m_name = name;
    m_length = length;

    m_first_base_in_kmer = -1;
    m_last_base_in_kmer = -1;

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
        if (length >= 16) {
            uint32_t kmer = kmers->starting_kmer_to_bits_forward(seq);
            for (int i = 15; i < length; ++i) {
                if (i > 15) {
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

    // See if the read failed any of the hard cut-offs.
    m_passed = true;
    if (args->min_length_set && m_length < args->min_length)
        m_passed = false;
    else if (args->max_length_set && m_length > args->max_length)
        m_passed = false;
    else if (args->min_mean_q_set && m_mean_quality < args->min_mean_q)
        m_passed = false;
    else if (args->min_window_q_set && m_window_quality < args->min_window_q)
        m_passed = false;

    m_first_base_in_kmer = -1;
    m_last_base_in_kmer = -1;
    if (!kmers->empty()) {
        for (int i = 0; i < length; ++i) {
            if (qualities[i] != 0) {
                if (m_first_base_in_kmer == -1)
                    m_first_base_in_kmer = i;
                m_last_base_in_kmer = i + 1;
            }
        }

        if (args->trim || args->split_set) {

            // Look at qualities to define 'bad ranges' of the read.
            if (args->split_set) {
                int i = 0;
                while (i < length) {
                    if (qualities[i] == 0.0) {
                        int bad_start = i;
                        while (i < length and qualities[i] == 0.0)
                            ++i;
                        int bad_end = i;
                        if (bad_end - bad_start >= args->split)
                            m_bad_ranges.push_back(std::pair<int,int>(bad_start, bad_end));
                    }
                    else
                        ++i;
                }
            }

            // If we're trimming, add the trimmed parts as bad ranges (unless they are already in there).
            if (args->trim) {
                if (m_first_base_in_kmer > 0) {
                    std::pair<int,int> trim_start(0, m_first_base_in_kmer);
                    if (m_bad_ranges.size() == 0 || m_bad_ranges.front() != trim_start)
                        m_bad_ranges.insert(m_bad_ranges.begin(), trim_start);
                }
                if (m_last_base_in_kmer != -1 && m_last_base_in_kmer < length) {
                    std::pair<int,int> trim_end(m_last_base_in_kmer, length);
                    if (m_bad_ranges.size() == 0 || m_bad_ranges.back() != trim_end)
                        m_bad_ranges.push_back(trim_end);
                }
            }

            if (m_bad_ranges.size() > 0) {
                int range_start = 0;
                int range_end;
                for (auto bad_range : m_bad_ranges) {
                    range_end = bad_range.first;
                    if (range_end - range_start > 0)
                        m_child_read_ranges.push_back(std::pair<int,int>(range_start, range_end));
                    range_start = bad_range.second;
                }
                range_end = length;
                if (range_end - range_start > 0)
                    m_child_read_ranges.push_back(std::pair<int,int>(range_start, range_end));
                for (size_t i = 0; i < m_child_read_ranges.size(); ++i) {
                    int child_start = m_child_read_ranges[i].first;
                    int child_length = m_child_read_ranges[i].second - child_start;
                    int child_end = child_start + child_length;
                    std::string child_name = m_name + "_" +
                            std::to_string(child_start+1) + "-" + std::to_string(child_end);
                    Read * child = new Read(child_name, seq + child_start, qscores + child_start, child_length,
                                            kmers, args);
                    m_child_reads.push_back(child);
                }
            }
        }
    }
}


Read::~Read() {
    for (auto child : m_child_reads)
        delete child;
}


std::string pad(std::string s, const size_t width)
{
    if (width > s.size())
        return s + std::string(width - s.size(), ' ');
    else
        return s;
}


std::string pad(int num, const size_t width)
{
    std::string s = std::to_string(num);
    return pad(s, width);
}


void Read::print_verbose_read_info() {
    std::cerr << "\n" << m_name << "\n";

    std::cerr << "            length = " << pad(m_length, 11);
    std::cerr << "mean quality = " << double_to_string(m_mean_quality);
    std::cerr << "      window quality = " << double_to_string(m_window_quality) << "\n";

    if (m_bad_ranges.size() > 0) {
        std::cerr << "        bad ranges = ";
        for (size_t i = 0; i < m_bad_ranges.size(); ++i) {
            std::cerr << m_bad_ranges[i].first << "-" << m_bad_ranges[i].second;
            if (i < m_bad_ranges.size() - 1)
                std::cerr << ", ";
        }
        std::cerr << "\n";
    }
    if (m_child_read_ranges.size() > 0) {
        std::cerr << "      child ranges = ";
        for (size_t i = 0; i < m_child_read_ranges.size(); ++i) {
            std::cerr << m_child_read_ranges[i].first << "-" << m_child_read_ranges[i].second;
            if (i < m_child_read_ranges.size() - 1)
                std::cerr << ", ";
        }
        std::cerr << "\n";
    }
    for (auto child : m_child_reads)
        child->print_verbose_read_info();
}


void Read::print_scores(size_t name_length) {
    std::cerr << pad(m_name, name_length) << "\t"
              << double_to_string(m_length_score) << "\t"
              << double_to_string(m_mean_quality) << "\t"
              << double_to_string(m_window_quality) << "\t"
              << double_to_string(m_final_score) << "\n";
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
        size_t i = j - window_size;
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
void Read::set_final_score(double length_weight, double mean_q_weight, double window_q_weight) {

    // First get the weighted geometric mean of the length score and the mean quality.
    double product = pow(m_length_score, length_weight) * pow(m_mean_quality, mean_q_weight);
    double total_weight = length_weight + mean_q_weight;
    double final_score = pow(product, 1.0 / total_weight);

    // Now scale that down using the ratio of window quality to mean quality.
    double scaling_factor;
    if (m_mean_quality > 0.0)
        scaling_factor = std::min(m_window_quality / m_mean_quality, 1.0);
    else
        scaling_factor = 1.0;
    total_weight = length_weight + mean_q_weight + window_q_weight;
    double window_weight_fraction = window_q_weight / total_weight;
    double non_window_weight_fraction = 1.0 - window_weight_fraction;
    scaling_factor = non_window_weight_fraction + (scaling_factor * window_weight_fraction);
    m_final_score = final_score * scaling_factor;
}


double Read::qscore_to_quality(char qscore) {
    int q = qscore - 33;
    return 1.0 - pow(10.0, -q / 10.0);
}
