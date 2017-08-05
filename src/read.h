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

#ifndef READ_H
#define READ_H


#include <string>
#include <unordered_set>
#include <vector>

#include "kmers.h"


class Read
{
public:
    Read(char * name, char * seq, char * qscores, int length, Kmers * kmers, int window_size,
         double length_weight, double mean_q_weight, double window_q_weight);

    void print_table_row();

    std::string m_name;

    int m_length;
    double m_mean_quality;
    double m_window_quality;

    double m_length_score;
    double m_mean_quality_score;
    double m_window_quality_score;
    double m_final_score;
    bool m_passed;

private:
    double get_mean_quality(std::vector<double> & qualities);
    double get_window_quality(std::vector<double> & qualities, size_t window_size);

    double get_length_score();
    double get_mean_quality_score();
    double get_window_quality_score();
    double get_final_score(double length_weight, double mean_q_weight, double window_q_weight);

    double qscore_to_quality(char qscore);
};


void print_read_table_header();

#endif // READ_H
