"""
Test unit suffix functionality for length-based options.

This module contains tests for the unit suffix parsing feature added to Filtlong.
To run them, execute `python3 -m unittest test.test_unit_suffixes` from the
root Filtlong directory.
"""

import unittest
import os
import subprocess


class TestUnitSuffixes(unittest.TestCase):

    def run_command(self, command):
        binary_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'bin', 'filtlong')
        input_path = os.path.join(os.path.dirname(__file__), 'test_sort.fastq')
        assembly_reference = os.path.join(os.path.dirname(__file__), 'test_reference.fasta')

        command = command.replace('filtlong', binary_path)
        command = command.replace('INPUT', input_path)
        command = command.replace('ASSEMBLY', assembly_reference)

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

    def test_target_bases_k_suffix(self):
        """Test that k suffix works for target_bases."""
        console_out, return_code = self.run_command('filtlong --target_bases 10k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 10' in console_out and 'bp' in console_out)

    def test_target_bases_kb_suffix(self):
        """Test that kb suffix works for target_bases."""
        console_out, return_code = self.run_command('filtlong --target_bases 10kb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 10' in console_out and 'bp' in console_out)

    def test_target_bases_g_suffix(self):
        """Test that g suffix works for target_bases."""
        console_out, return_code = self.run_command('filtlong --target_bases 1g INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 1' in console_out and '000' in console_out and 'bp' in console_out)

    def test_target_bases_gb_suffix(self):
        """Test that gb suffix works for target_bases."""
        console_out, return_code = self.run_command('filtlong --target_bases 1gb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 1' in console_out and '000' in console_out and 'bp' in console_out)

    def test_target_bases_m_suffix(self):
        """Test that m suffix works for target_bases."""
        console_out, return_code = self.run_command('filtlong --target_bases 0.01m INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        self.assertTrue('target: 10,000 bp' in console_out)

    def test_target_bases_mb_suffix(self):
        """Test that mb suffix works for target_bases."""
        console_out, return_code = self.run_command('filtlong --target_bases 0.01mb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 10' in console_out and 'bp' in console_out)

    def test_target_bases_case_insensitive(self):
        """Test that suffixes are case insensitive."""
        console_out, return_code = self.run_command('filtlong --target_bases 10K INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 10' in console_out and 'bp' in console_out)

        console_out, return_code = self.run_command('filtlong --target_bases 10KB INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 10' in console_out and 'bp' in console_out)

        console_out, return_code = self.run_command('filtlong --target_bases 0.01M INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 10' in console_out and 'bp' in console_out)

        console_out, return_code = self.run_command('filtlong --target_bases 0.01MB INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 10' in console_out and 'bp' in console_out)

    def test_problem_statement_examples(self):
        """Test the specific examples from the problem statement."""
        # Test 3.5mb should become 3500000
        console_out, return_code = self.run_command('filtlong --target_bases 3.5mb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 3' in console_out and '500' in console_out and 'bp' in console_out)

        # Test 1kb should become 1000
        console_out, return_code = self.run_command('filtlong --target_bases 1kb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 1' in console_out and '000' in console_out and 'bp' in console_out)

        # Test 1k should become 1000
        console_out, return_code = self.run_command('filtlong --target_bases 1k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Check for target value (may have commas depending on locale)
        self.assertTrue('target: 1' in console_out and '000' in console_out and 'bp' in console_out)

    def test_min_length_suffix(self):
        """Test that unit suffixes work for min_length."""
        console_out, return_code = self.run_command('filtlong --min_length 1k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Should filter out reads shorter than 1000bp (all test reads are 5000bp)
        # so all should pass through

    def test_min_length_g_suffix(self):
        """Test that g suffix works for min_length."""
        console_out, return_code = self.run_command('filtlong --min_length 1g INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Should filter out all reads since they're much shorter than 1GB

    def test_max_length_suffix(self):
        """Test that unit suffixes work for max_length."""
        console_out, return_code = self.run_command('filtlong --max_length 10k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Should pass all reads since they're all 5000bp

    def test_max_length_gb_suffix(self):
        """Test that gb suffix works for max_length."""
        console_out, return_code = self.run_command('filtlong --max_length 1gb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Should pass all reads since they're much shorter than 1GB

    def test_split_suffix(self):
        """Test that unit suffixes work for split."""
        console_out, return_code = self.run_command('filtlong -a ASSEMBLY --split 1k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)

    def test_split_g_suffix(self):
        """Test that g suffix works for split."""
        console_out, return_code = self.run_command('filtlong -a ASSEMBLY --split 1g INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)

    def test_invalid_suffix_error(self):
        """Test that invalid suffixes produce appropriate error messages."""
        console_out, return_code = self.run_command('filtlong --target_bases 10xyz INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 1)
        self.assertTrue('invalid value' in console_out)

    def test_no_numeric_value_error(self):
        """Test that suffixes without numeric values produce appropriate error messages."""
        console_out, return_code = self.run_command('filtlong --target_bases k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 1)
        self.assertTrue('invalid value' in console_out)

    def test_decimal_values_work(self):
        """Test that decimal values work with suffixes."""
        console_out, return_code = self.run_command('filtlong --target_bases 3.5k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        self.assertTrue('target: 3,500 bp' in console_out)

    def test_backwards_compatibility(self):
        """Test that plain numbers still work (backwards compatibility)."""
        console_out, return_code = self.run_command('filtlong --target_bases 10000 INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        self.assertTrue('target: 10,000 bp' in console_out)

    def test_negative_values_with_suffixes(self):
        """Test that negative values with suffixes are properly rejected."""
        console_out, return_code = self.run_command('filtlong --target_bases -10k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 1)
        self.assertTrue('Error: the value for --target_bases must be a positive integer' in console_out)

        console_out, return_code = self.run_command('filtlong --min_length -5kb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 1)
        self.assertTrue('Error: the value for --min_length must be a positive integer' in console_out)

    def test_min_length_short_option_suffix(self):
        """Test that unit suffixes work for -l option."""
        console_out, return_code = self.run_command('filtlong -l 5k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Should pass all reads since they're all 5000bp

    def test_max_length_short_option_suffix(self):
        """Test that unit suffixes work for -L option."""
        console_out, return_code = self.run_command('filtlong -L 10k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 0)
        # Should pass all reads since they're all 5000bp

    def test_negative_values_with_suffixes_short_options(self):
        """Test that negative values with suffixes are properly rejected for short options."""
        console_out, return_code = self.run_command('filtlong -l -10k INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 1)
        self.assertTrue('Error: the value for --min_length must be a positive integer' in console_out)

        console_out, return_code = self.run_command('filtlong -L -5kb INPUT > OUTPUT.fastq')
        self.assertEqual(return_code, 1)
        self.assertTrue('Error: the value for --max_length must be a positive integer' in console_out)
