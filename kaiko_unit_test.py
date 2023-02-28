import unittest
import unit_test_util as util
import glob
import re

import pandas as pd

from pathlib import Path

unit_test_dir = "./Kaiko_volume/Kaiko_stationary_files/unit_test_correct_output/"


class TestStringMethods(unittest.TestCase):

    def test_params_1(self):
        topk_test, beam_size_test = False, 5
        dir_ending = "mgf_large_unit_test_topk_{a1}_beam_size_{a3}".format(a1 = topk_test, a3 = beam_size_test)

        test_output = Path("./Kaiko_volume/Kaiko_intermediate/denovo_output/" + dir_ending + "/*unit_test_prop_*.txt")
        correct_output = Path(unit_test_dir + dir_ending + "/*unit_test_prop_*.txt") 
        correct_output = glob.glob(str(correct_output.absolute()))
        print("Comparing to {a1} unit test files. \n".format(a1 = len(correct_output)))
        print("Running test.\n")

        util.run_new_parameters(topk = topk_test, beam_size = beam_size_test)

        test_output = glob.glob(str(test_output.absolute()))

        for i, denovo_out in enumerate(test_output):
            unit_output = Path(correct_output[i])
            denovo_out = Path(denovo_out)
            print(denovo_out.name)
            self.assertEqual(denovo_out.name, unit_output.name)

            with denovo_out.open() as denovo_f, unit_output.open() as unit_f:
                for new_line, unit_line in zip(denovo_f, unit_f):
                    self.assertEqual(new_line, unit_line)

    def test_params_2(self):
        topk_test, beam_size_test = False, 7
        dir_ending = "mgf_large_unit_test_topk_{a1}_beam_size_{a3}".format(a1 = topk_test, a3 = beam_size_test)

        test_output = Path("./Kaiko_volume/Kaiko_intermediate/denovo_output/" + dir_ending + "/*unit_test_prop_*.txt")
        correct_output = Path(unit_test_dir + dir_ending + "/*unit_test_prop_*.txt") 
        correct_output = glob.glob(str(correct_output.absolute()))
        print("Comparing to {a1} unit test files. \n".format(a1 = len(correct_output)))
        print("Running test.\n")

        util.run_new_parameters(topk = topk_test, beam_size = beam_size_test)

        test_output = glob.glob(str(test_output.absolute()))

        for i, denovo_out in enumerate(test_output):
            unit_output = Path(correct_output[i])
            denovo_out = Path(denovo_out)
            print(denovo_out.name)
            self.assertEqual(denovo_out.name, unit_output.name)

            with denovo_out.open() as denovo_f, unit_output.open() as unit_f:
                for new_line, unit_line in zip(denovo_f, unit_f):
                    self.assertEqual(new_line, unit_line)

    
    def test_params_topk_1(self):
        topk_test, beam_size_test = True, 5
        dir_ending = "mgf_large_unit_test_topk_{a1}_beam_size_{a3}".format(a1 = topk_test, a3 = beam_size_test)

        test_output = Path("./Kaiko_volume/Kaiko_intermediate/denovo_output/" + dir_ending + "/*unit_test_prop_*.txt")
        correct_output = Path(unit_test_dir + dir_ending + "/*unit_test_prop_*.txt") 
        correct_output = glob.glob(str(correct_output))
        print("Comparing to {a1} unit test files. \n".format(a1 = len(correct_output)))
        print("Running test.\n")

        util.run_new_parameters(topk = topk_test, beam_size = beam_size_test)

        test_output = glob.glob(str(test_output))

        for i, denovo_out in enumerate(test_output):
            unit_output = Path(correct_output[i])
            denovo_out = Path(denovo_out)
            print(denovo_out.name)
            self.assertEqual(denovo_out.name, unit_output.name)

            with denovo_out.open() as denovo_f, unit_output.open() as unit_f:
                old_scan = ""
                new_lines = []
                unit_lines = []
                for new_line, unit_line in zip(denovo_f, unit_f):
                    scan = re.sub("(^[0-9]+\\:[0-9]+).*$", "\\1", new_line)
                    if old_scan == scan:
                        new_lines += [new_line]
                        unit_lines += [unit_line]
                    else:
                        old_scan = scan
                        self.assertEqual(sorted(new_lines), sorted(unit_lines))
                        new_lines = [new_line]
                        unit_lines = [unit_line]


    def test_params_topk_2(self):
        topk_test, beam_size_test = True, 7
        dir_ending = "mgf_large_unit_test_topk_{a1}_beam_size_{a3}".format(a1 = topk_test, a3 = beam_size_test)

        test_output = Path("./Kaiko_volume/Kaiko_intermediate/denovo_output/" + dir_ending + "/*unit_test_prop_*.txt")
        correct_output = Path(unit_test_dir + dir_ending + "/*unit_test_prop_*.txt") 
        correct_output = glob.glob(str(correct_output))
        print("Comparing to {a1} unit test files. \n".format(a1 = len(correct_output)))
        print("Running test.\n")

        util.run_new_parameters(topk = topk_test, beam_size = beam_size_test)

        test_output = glob.glob(str(test_output))

        for i, denovo_out in enumerate(test_output):
            unit_output = Path(correct_output[i])
            denovo_out = Path(denovo_out)
            print(denovo_out.name)
            self.assertEqual(denovo_out.name, unit_output.name)

            with denovo_out.open() as denovo_f, unit_output.open() as unit_f:
                old_scan = ""
                new_lines = []
                unit_lines = []
                for new_line, unit_line in zip(denovo_f, unit_f):
                    scan = re.sub("(^[0-9]+\\:[0-9]+).*$", "\\1", new_line)
                    if old_scan == scan:
                        new_lines += [new_line]
                        unit_lines += [unit_line]
                    else:
                        old_scan = scan
                        self.assertEqual(sorted(new_lines), sorted(unit_lines))
                        new_lines = [new_line]
                        unit_lines = [unit_line]

        

if __name__ == '__main__':
    unittest.main()










