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

#ifndef READ_H
#define READ_H


#include <string>
#include <unordered_set>
#include <utility>
#include <tuple>

#include "kmers.h"
#include "arguments.h"


class Read
{
public:
    Read(std::string name, char * seq, char * qscores, int length, Kmers * kmers, Arguments * args);
    ~Read();

    void print_verbose_read_info();
    void print_scores(size_t name_length);

    void set_final_score(double length_weight, double mean_q_weight, double window_q_weight);

    std::string m_name;

    int m_length;
    double m_length_score;

    double m_mean_quality;
    double m_window_quality;

    double m_final_score;
    bool m_passed;

    int m_first_base_in_kmer;
    int m_last_base_in_kmer;
    std::vector<std::pair<int,int> > m_bad_ranges;

    std::vector<Read *> m_child_reads;
    std::vector<std::pair<int,int> > m_child_read_ranges;

private:
    double get_mean_quality(std::vector<double> & qualities);
    double get_window_quality(std::vector<double> & qualities, size_t window_size);

    double get_length_score();

    double qscore_to_quality(char qscore);
};


#endif // READ_H
