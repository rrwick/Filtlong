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


#include "arguments.h"

#include <iostream>
#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>

#include "args.h"


struct DoublesReader
{
    void operator()(const std::string &name, const std::string &value, double &destination) {
        try {
            if (value.find_first_not_of("0123456789.") != std::string::npos)
                throw std::invalid_argument("");
            destination = std::stod(value);
        }
        catch ( ... ) {
            std::ostringstream problem;
            problem << "Error: argument '" << name << "' received invalid value type '" << value << "'";
            throw args::ParseError(problem.str());
        }
    }
};


typedef args::ValueFlag<double, DoublesReader> d_arg;
typedef args::ValueFlag<long long> i_arg;
typedef args::ValueFlag<std::string> s_arg;
typedef args::Flag f_arg;


Arguments::Arguments(int argc, char **argv) {

    args::ArgumentParser parser("Filtlong: a quality filtering tool for Nanopore and PacBio reads",
                                "For more information, go to: https://github.com/rrwick/Filtlong");
    parser.LongSeparator(" ");

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    int terminal_width = w.ws_col;

    int indent_size;
    if (terminal_width > 120)
        indent_size = 4;
    else if (terminal_width > 80)
        indent_size = 3;
    else if (terminal_width > 60)
        indent_size = 2;
    else
        indent_size = 1;

    parser.helpParams.showTerminator = false;
    parser.helpParams.progindent = 0;
    parser.helpParams.descriptionindent = 0;
    parser.helpParams.width = terminal_width;
    parser.helpParams.flagindent = indent_size;
    parser.helpParams.eachgroupindent = indent_size;

    args::Positional<std::string> input_reads_arg(parser, "input_reads",
                                      "input long reads to be filtered");

    args::Group thresholds_group(parser, "output thresholds:");
    i_arg target_bases_arg(thresholds_group, "int",
                           "keep only the best reads up to this many total bases",
                           {'t', "target_bases"});
    d_arg keep_percent_arg(thresholds_group, "float",
                           "keep only this percentage of the best reads (measured by bases)",
                           {'p', "keep_percent"});
    i_arg min_length_arg(thresholds_group, "int",
                         "minimum length threshold",
                         {"min_length"});
    i_arg max_length_arg(thresholds_group, "int",
                         "maximum length threshold",
                         {"max_length"});
    d_arg min_mean_q_arg(thresholds_group, "float",
                         "minimum mean quality threshold",
                         {"min_mean_q"});
    d_arg min_window_q_arg(thresholds_group, "float",
                           "minimum window quality threshold",
                           {"min_window_q"});

    args::Group references_group(parser, "NLexternal references "   // The NL at the start results in a newline
            "(if provided, read quality will be determined using these instead of from the Phred scores):");
    s_arg assembly_arg(references_group, "file",
                       "reference assembly in FASTA format",
                        {'a', "assembly"});
    s_arg illumina_1_arg(references_group, "file",
                         "reference Illumina reads in FASTQ format",
                         {'1', "illumina_1"});
    s_arg illumina_2_arg(references_group, "file",
                         "reference Illumina reads in FASTQ format",
                         {'2', "illumina_2"});

    args::Group score_weights_group(parser, "NLscore weights "    // The NL at the start results in a newline
                                            "(control the relative contribution of each score to the final read score):");
    d_arg length_weight_arg(score_weights_group, "float",
                            "weight given to the length score (default: 1)",
                            {"length_weight"}, 1.0);
    d_arg mean_q_weight_arg(score_weights_group, "float",
                            "weight given to the mean quality score (default: 1)",
                            {"mean_q_weight"}, 1.0);
    d_arg window_q_weight_arg(score_weights_group, "float",
                              "weight given to the window quality score (default: 1)",
                              {"window_q_weight"}, 1.0);

    // TO DO: half_length_score: reads of this length get a 50% score (currently hard-coded at 5000)

    args::Group manipulation_group(parser, "NLread manipulation:");    // The NL at the start results in a newline
    f_arg trim_arg(manipulation_group, "trim",
                   "trim non-k-mer-matching bases from start/end of reads",
                   {"trim"});
    i_arg split_arg(manipulation_group, "split",
                    "split reads at this many (or more) consecutive non-k-mer-matching bases",
                    {"split"});

    args::Group other_group(parser, "NLother:");    // The NL at the start results in a newline
    i_arg window_size_arg(other_group, "int",
                          "size of sliding window used when measuring window quality (default: 250)",
                          {"window_size"}, 250);
    f_arg verbose_arg(other_group, "verbose",
                      "verbose output to stderr with info for each read",
                      {"verbose"});
    f_arg version_arg(other_group, "version",
                      "display the program version and quit",
                      {"version"});

    args::HelpFlag help(parser, "help",
                        "display this help menu",
                        {'h', "help"});


    parsing_result = GOOD;
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        std::cerr << parser;
        parsing_result = HELP;
        return;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << "\n";
        parsing_result = BAD;
        return;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << "\n";
        parsing_result = BAD;
        return;
    }
    if (argc == 1) {
        std::cerr << parser;
        parsing_result = HELP;
        return;
    }
    if (args::get(version_arg)) {
        parsing_result = VERSION;
        return;
    }

    input_reads = args::get(input_reads_arg);
    if (input_reads.empty()) {
        std::cerr << "Error: input reads are required" << "\n";
        parsing_result = BAD;
        return;
    }

    target_bases_set = bool(target_bases_arg);
    target_bases = args::get(target_bases_arg);

    keep_percent_set = bool(keep_percent_arg);
    keep_percent = args::get(keep_percent_arg);

    assembly_set = bool(assembly_arg);
    assembly = args::get(assembly_arg);

    if (bool(illumina_1_arg))
        illumina_reads.push_back(args::get(illumina_1_arg));
    if (bool(illumina_2_arg))
        illumina_reads.push_back(args::get(illumina_2_arg));

    min_length_set = bool(min_length_arg);
    min_length = args::get(min_length_arg);

    max_length_set = bool(max_length_arg);
    max_length = args::get(max_length_arg);

    min_mean_q_set = bool(min_mean_q_arg);
    min_mean_q = args::get(min_mean_q_arg);

    min_window_q_set = bool(min_window_q_arg);
    min_window_q = args::get(min_window_q_arg);

    length_weight = args::get(length_weight_arg);
    mean_q_weight = args::get(mean_q_weight_arg);
    window_q_weight = args::get(window_q_weight_arg);

    trim = args::get(trim_arg);

    split_set = bool(split_arg);
    split = args::get(split_arg);

    window_size = args::get(window_size_arg);
    verbose = args::get(verbose_arg);

    bool some_reference = (illumina_reads.size() > 0 || assembly_set);
    if (trim && !some_reference) {
        std::cerr << "Error: assembly or read reference is required to use --trim" << "\n";
        parsing_result = BAD;
        return;
    }
    if (split_set && !some_reference) {
        std::cerr << "Error: assembly or read reference is required to use --split" << "\n";
        parsing_result = BAD;
        return;
    }

    // Check to make sure files exist.
    std::vector<std::string> files;
    files.push_back(input_reads);
    for (auto f : illumina_reads)
        files.push_back(f);
    if (assembly_set)
        files.push_back(assembly);
    for (auto f : files) {
        if (!does_file_exist(f)) {
            std::cerr << "Error: cannot find file: " << f << "\n";
            parsing_result = BAD;
            return;
        }
    }

    // If nothing is set, then Filtlong won't do anything. Give an error message and quit.
    if (!trim && !split_set && !target_bases_set && !keep_percent_set &&
            !min_length_set && !max_length_set && !min_mean_q_set && !min_window_q_set) {
        std::cerr << "Error: no thresholds set, you must use one of the following options:\n";
        std::cerr << "target_bases, keep_percent, min_length, max_length, min_mean_q, min_window_q, trim, split\n";
        parsing_result = BAD;
        return;
    }

    // Non-positive target_bases doesn't make sense.
    if (target_bases_set && target_bases <= 0) {
        std::cerr << "Error: the value for --target_bases must be a positive integer\n";
        parsing_result = BAD;
        return;
    }

    // Non-positive min_length doesn't make sense.
    if (min_length_set && min_length <= 0) {
        std::cerr << "Error: the value for --min_length must be a positive integer\n";
        parsing_result = BAD;
        return;
    }
 
    // Non-positive max_length doesn't make sense.
    if (max_length_set && max_length <= 0) {
        std::cerr << "Error: the value for --max_length must be a positive integer\n";
        parsing_result = BAD;
        return;
    }

    // keep_percent must be between 0 and 100 (exclusive).
    if (keep_percent_set && (keep_percent <= 0.0 || keep_percent >= 100.0)) {
        std::cerr << "Error: the value for --keep_percent must be greater than 0 and less than 100\n";
        parsing_result = BAD;
        return;
    }

    // min_mean_q and min_window_q must be positive.
    if (min_mean_q_set && min_mean_q <= 0.0) {
        std::cerr << "Error: the value for --min_mean_q must be greater than 0\n";
        parsing_result = BAD;
        return;
    }
    if (min_window_q_set && min_window_q <= 0.0) {
        std::cerr << "Error: the value for --min_window_q must be greater than 0\n";
        parsing_result = BAD;
        return;
    }

    // Negative weights don't make sense.
    if (length_weight < 0.0 || mean_q_weight < 0.0 || window_q_weight < 0.0) {
        std::cerr << "Error: weight values cannot be negative\n";
        parsing_result = BAD;
        return;
    }

    // Non-positive split doesn't make sense.
    if (split_set && split <= 0) {
        std::cerr << "Error: the value for --split must be a positive integer\n";
        parsing_result = BAD;
        return;
    }

    // Non-positive window_size doesn't make sense.
    if (window_size <= 0) {
        std::cerr << "Error: the value for --window_size must be a positive integer\n";
        parsing_result = BAD;
        return;
    }
}


bool Arguments::does_file_exist(std::string filename){
    std::ifstream infile(filename);
    return infile.good();
}
