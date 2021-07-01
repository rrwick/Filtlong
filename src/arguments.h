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

    bool target_bases_set;
    long long target_bases;

    bool keep_percent_set;
    double keep_percent;

    bool min_length_set;
    int min_length;

    bool max_length_set;
    int max_length;

    bool min_mean_q_set;
    double min_mean_q;

    bool min_window_q_set;
    double min_window_q;

    bool assembly_set;
    std::string assembly;
    std::vector<std::string> illumina_reads;

    double length_weight;
    double mean_q_weight;
    double window_q_weight;

    bool trim;

    bool split_set;
    int split;

    int window_size;
    bool verbose;


private:
    bool does_file_exist(std::string fileName);
};

#endif // ARGUMENTS_H
