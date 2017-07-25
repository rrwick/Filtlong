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

#include "kmers.h"


class Read
{
public:
    Read(std::string name, std::string seq, std::string quals, Kmers * kmers, int window_size);

    std::string m_name;
    int m_length;
    double m_mean_quality;
    double m_window_quality;

private:
    double get_mean_quality(std::string * quals);
    double get_window_quality(std::string * quals, int window_size);
};

#endif // READ_H
