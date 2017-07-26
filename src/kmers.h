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

#ifndef KMERS_H
#define KMERS_H


#include <string>
#include <vector>
#include <unordered_set>


class Kmers
{
public:
    Kmers();

    void add_read_fastqs(std::vector<std::string> filenames);
    void add_assembly_fasta(std::string filename);
    bool is_kmer_present(uint32_t kmer);

private:
    std::unordered_set<uint32_t> m_kmers;
    std::unordered_set<uint32_t> m_first_time_kmers;

    int add_reference(std::string filename, bool require_two_kmer_copies);
    void add_kmer_require_one_copy(uint32_t kmer);
    void add_kmer_require_two_copies(uint32_t kmer);

    uint32_t starting_kmer_to_bits_forward(char * sequence);
    uint32_t starting_kmer_to_bits_reverse(char * sequence);

    uint32_t base_to_bits_forward(char base);
    uint32_t base_to_bits_reverse(char base);
};


#endif // KMERS_H
