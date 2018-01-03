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


class TestErrorMessages(unittest.TestCase):

    def run_command(self, command):
        binary_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'bin', 'filtlong')
        input_path = os.path.join(os.path.dirname(__file__), 'test_sort.fastq')
        input_fasta = os.path.join(os.path.dirname(__file__), 'test_sort.fasta')
        bad_fastq = os.path.join(os.path.dirname(__file__), 'test_bad_fastq.fastq')
        assembly_reference = os.path.join(os.path.dirname(__file__), 'test_reference.fasta')
        illumina_reference_1 = os.path.join(os.path.dirname(__file__), 'test_reference_1.fastq.gz')
        illumina_reference_2 = os.path.join(os.path.dirname(__file__), 'test_reference_2.fastq.gz')

        command = command.replace('filtlong', binary_path)
        command = command.replace('INPUT', input_path)
        command = command.replace('FASTA', input_fasta)
        command = command.replace('BADFASTQ', bad_fastq)
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
        return err.decode(), p.returncode

    def tearDown(self):
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

    def test_empty_command(self):
        """
        If no arguments are given at all, it should display the help text.
        """
        console_out, return_code = self.run_command('filtlong')
        self.assertTrue('usage:' in console_out)
        self.assertTrue('Filtlong:' in console_out)
        self.assertEqual(return_code, 0)

    def test_no_thresholds(self):
        console_out, return_code = self.run_command('filtlong INPUT > OUTPUT.fastq')
        self.assertTrue('Error: no thresholds set' in console_out)
        self.assertEqual(return_code, 1)

    def test_input_reads_do_not_exit(self):
        console_out, return_code = self.run_command('filtlong --target_bases 1000 BAD_FILENAME > OUTPUT.fastq')
        self.assertTrue('Error: cannot find file' in console_out)
        self.assertEqual(return_code, 1)

    def test_reference_assembly_does_not_exit(self):
        console_out, return_code = self.run_command('filtlong -a BAD_FILENAME --target_bases 1000 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: cannot find file' in console_out)
        self.assertEqual(return_code, 1)

    def test_reference_reads_1_does_not_exit(self):
        console_out, return_code = self.run_command('filtlong -1 BAD_FILENAME -2 ILLUMINA_2 --target_bases 5000 '
                                                    'INPUT > OUTPUT.fastq')
        self.assertTrue('Error: cannot find file' in console_out)
        self.assertEqual(return_code, 1)

    def test_reference_reads_2_does_not_exit(self):
        console_out, return_code = self.run_command('filtlong -1 ILLUMINA_1 -2 BAD_FILENAME --target_bases 5000 '
                                                    'INPUT > OUTPUT.fastq')
        self.assertTrue('Error: cannot find file' in console_out)
        self.assertEqual(return_code, 1)

    def test_target_bases_too_low_1(self):
        console_out, return_code = self.run_command('filtlong --target_bases 0 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --target_bases must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_target_bases_too_low_2(self):
        console_out, return_code = self.run_command('filtlong --target_bases -10 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --target_bases must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_keep_percent_too_low(self):
        console_out, return_code = self.run_command('filtlong --keep_percent 0 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --keep_percent must be greater than 0 and less than 100' in console_out)
        self.assertEqual(return_code, 1)

    def test_keep_percent_too_high_1(self):
        console_out, return_code = self.run_command('filtlong --keep_percent 100 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --keep_percent must be greater than 0 and less than 100' in console_out)
        self.assertEqual(return_code, 1)

    def test_keep_percent_too_high_2(self):
        console_out, return_code = self.run_command('filtlong --keep_percent 111.1 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --keep_percent must be greater than 0 and less than 100' in console_out)
        self.assertEqual(return_code, 1)

    def test_min_length_too_low_1(self):
        console_out, return_code = self.run_command('filtlong --min_length 0 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --min_length must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_min_length_too_low_2(self):
        console_out, return_code = self.run_command('filtlong --min_length -10 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --min_length must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_min_mean_q_too_low(self):
        console_out, return_code = self.run_command('filtlong --min_mean_q 0 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --min_mean_q must be greater than 0' in console_out)
        self.assertEqual(return_code, 1)

    def test_min_window_q_too_low(self):
        console_out, return_code = self.run_command('filtlong --min_window_q 0 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --min_window_q must be greater than 0' in console_out)
        self.assertEqual(return_code, 1)

    def test_trim_no_reference(self):
        console_out, return_code = self.run_command('filtlong --trim INPUT > OUTPUT.fastq')
        self.assertTrue('Error: assembly or read reference is required to use --trim' in console_out)
        self.assertEqual(return_code, 1)

    def test_split_no_reference(self):
        console_out, return_code = self.run_command('filtlong --split 250 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: assembly or read reference is required to use --split' in console_out)
        self.assertEqual(return_code, 1)

    def test_split_too_low_1(self):
        console_out, return_code = self.run_command('filtlong -a ASSEMBLY --split 0 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --split must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_split_too_low_2(self):
        console_out, return_code = self.run_command('filtlong -a ASSEMBLY --split -10 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --split must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_window_size_too_low_1(self):
        console_out, return_code = self.run_command('filtlong --min_length 1000 --window_size 0 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --window_size must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_window_size_too_low_2(self):
        console_out, return_code = self.run_command('filtlong --min_length 1000 --window_size -10 INPUT > OUTPUT.fastq')
        self.assertTrue('Error: the value for --window_size must be a positive integer' in console_out)
        self.assertEqual(return_code, 1)

    def test_fasta_input(self):
        console_out, return_code = self.run_command('filtlong --target_bases 1000 FASTA > OUTPUT.fastq')
        self.assertTrue('Error: FASTA input not supported without an external reference' in console_out)
        self.assertEqual(return_code, 1)

    def test_bad_fastq(self):
        console_out, return_code = self.run_command('filtlong --target_bases 1000 BADFASTQ > OUTPUT.fastq')
        self.assertTrue('Error: incorrect FASTQ format for read' in console_out)
        self.assertEqual(return_code, 1)
