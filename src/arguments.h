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

#ifndef ARGUMENTS_H
#define ARGUMENTS_H


#include <string>
#include <vector>


enum ParsingResult {GOOD, BAD, HELP, VERSION};


class Arguments
{
public:
    Arguments(int argc, char **argv);

    ParsingResult parsing_result;

    std::string input_reads;

    double min_score;
    long long target_bases;
    double keep_percent;
    std::string assembly;
    std::vector<std::string> illumina_reads;
    int min_length;
    double min_mean_q;
    double min_window_q;
    double length_weight;
    double mean_q_weight;
    double window_q_weight;
    int window_size;
    bool verbose;

    bool min_score_set;
    bool target_bases_set;
    bool keep_percent_set;
    bool assembly_set;
    bool min_length_set;
    bool min_mean_q_set;
    bool min_window_q_set;
};

#endif // ARGUMENTS_H
