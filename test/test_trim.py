"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Filtlong

This module contains some tests for Filtlong. To run them, execute `python3 -m unittest` from the
root Filtlong directory.

This file is part of Filtlong. Filtlong is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Filtlong is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Filtlong. If
not, see <http://www.gnu.org/licenses/>.
"""

import unittest
import os
import subprocess


def load_fastq(filename):
    reads = []
    with open(filename, 'rb') as fastq:
        for line in fastq:
            stripped_line = line.strip()
            if len(stripped_line) == 0:
                continue
            if not stripped_line.startswith(b'@'):
                continue
            name = stripped_line[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            reads.append((name, sequence, qualities))
    return reads


class TestTrim(unittest.TestCase):

    def run_command(self, command):
        binary_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'bin', 'filtlong')
        input_path = os.path.join(os.path.dirname(__file__), 'test_trim.fastq')
        assembly_reference = os.path.join(os.path.dirname(__file__), 'test_reference.fasta')
        illumina_reference_1 = os.path.join(os.path.dirname(__file__), 'test_reference_1.fastq.gz')
        illumina_reference_2 = os.path.join(os.path.dirname(__file__), 'test_reference_2.fastq.gz')

        command = command.replace('filtlong', binary_path)
        command = command.replace('INPUT', input_path)
        command = command.replace('ASSEMBLY', assembly_reference)
        command = command.replace('ILLUMINA_1', illumina_reference_1)
        command = command.replace('ILLUMINA_2', illumina_reference_2)

        output_name = 'TEMP_' + str(os.getpid())
        command = command.replace('OUTPUT', output_name)
        try:
            self.output_file = [x for x in command.split() if output_name in x][0]
        except IndexError:
            self.output_file = ''
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        _, err = p.communicate()
        return err.decode()

    def tearDown(self):
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

    def check_one_read(self, read, length, start_seq, end_seq):
        self.assertEqual(len(read[1]), length)
        self.assertEqual(len(read[2]), length)
        self.assertTrue(read[1].startswith(start_seq))
        self.assertTrue(read[1].endswith(end_seq))

    def test_trim_1(self):
        """
        Tests read trimming console output using an assembly reference.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --trim INPUT > OUTPUT.fastq')
        self.assertTrue('4 reads (4,901 bp)' in console_out)
        self.assertTrue('after trimming: 4 reads (4,824 bp)' in console_out)

    def test_trim_2(self):
        """
        Tests read trimming console output using Illumina read references.
        """
        console_out = self.run_command('filtlong -1 ILLUMINA_1 -2 ILLUMINA_2 --trim INPUT > OUTPUT.fastq')
        self.assertTrue('4 reads (4,901 bp)' in console_out)
        self.assertTrue('after trimming: 4 reads (4,824 bp)' in console_out)

    def test_trim_3(self):
        """
        Actually loads the trimmed reads to verify that they look good.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --trim INPUT > OUTPUT.fastq')
        trimmed_reads = load_fastq(self.output_file)
        self.assertEqual(len(trimmed_reads), 4)
        self.check_one_read(trimmed_reads[0], 1300, b'GCCCTGGC', b'GGGTCCAG')
        self.check_one_read(trimmed_reads[1], 701 - 20, b'GATTTATA', b'ATGGCGAC')
        self.check_one_read(trimmed_reads[2], 1000 - 30, b'CTTGAACA', b'TCCTCCAG')
        self.check_one_read(trimmed_reads[3], 1900 - 27, b'CCTTTCTT', b'TGATCACC')

    def test_trim_names(self):
        """
        Makes sure the trimmed reads are correctly named.
        """
        console_out = self.run_command('filtlong -a ASSEMBLY --trim INPUT > OUTPUT.fastq')
        trimmed_reads = load_fastq(self.output_file)
        self.assertEqual(len(trimmed_reads), 4)
        self.assertEqual(trimmed_reads[0][0], b'test_trim_1')
        self.assertEqual(trimmed_reads[1][0], b'test_trim_2_21-701')
        self.assertEqual(trimmed_reads[2][0], b'test_trim_3_1-970')
        self.assertEqual(trimmed_reads[3][0], b'test_trim_4_13-1885')
