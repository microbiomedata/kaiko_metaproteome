import unittest
import unit_test_util as util
import glob
import re

from pathlib import Path

unit_test_dir = "./Kaiko_volume/Kaiko_stationary_files/unit_test_correct_output/"
score_cutoff = 3.00

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
            self.assertEqual(denovo_out.name, unit_output.name)

            with denovo_out.open() as denovo_f, unit_output.open() as unit_f:
                print(denovo_out.name)
                next(denovo_f), next(unit_f)
                old_scan = ""
                new_lines = dict()
                unit_lines = dict()
                for new_line in denovo_f:
                    new_scan = re.sub("(^[0-9]+\\:[0-9]+).*$", "\\1", new_line)
                    new_score = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\2", new_line)
                    if new_score == 'f\n':
                        new_score = 'inf'
                    new_score = float(new_score)
                    new_line = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\1\\3", new_line)

                    if old_scan == new_scan:
                        if new_line not in new_lines.keys():
                            new_lines[new_line] = new_score
                        elif new_lines[new_line] > new_score:
                            new_lines[new_line] = new_score
                    else:
                        unit_scan = old_scan
                        while (unit_scan == old_scan):
                            unit_line = unit_f.readline()
                            unit_scan = re.sub("(^[0-9]+\\:[0-9]+).*$", "\\1", unit_line)
                            unit_score = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\2", unit_line)
                            if unit_score == 'f\n':
                                unit_score = 'inf'
                            unit_score = float(unit_score)
                            unit_line = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\1\\3", unit_line)

                            if unit_scan == old_scan:
                                if unit_line not in unit_lines.keys():
                                    unit_lines[unit_line] = unit_score
                                elif unit_lines[unit_line] > unit_score:
                                    unit_lines[unit_line] = unit_score
                        old_scan = new_scan

                        new_lines_conf = {k: new_lines[k] for k in new_lines if new_lines[k] < score_cutoff}
                        unit_lines_conf = {k: unit_lines[k] for k in unit_lines if unit_lines[k] < score_cutoff}

                        ## Checking that predictions with a good enough score are in both results.
                        self.assertTrue(set(new_lines_conf.keys()).issubset(unit_lines.keys()))
                        self.assertTrue(set(unit_lines_conf.keys()).issubset(new_lines.keys()))

                        ## Checking that score of good predictions is within no more than 0.05 when compared to unit test.
                        for k in unit_lines_conf.keys():
                            self.assertAlmostEqual((unit_lines_conf[k] - new_lines[k]) * 2, 0, places = 1)
                        for k in new_lines_conf.keys():
                            self.assertAlmostEqual((unit_lines[k] - new_lines_conf[k]) * 2, 0, places = 1)
                        new_lines = dict()
                        unit_lines = dict()
                        new_lines[new_line] = new_score
                        unit_lines[unit_line] = unit_score


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
            self.assertEqual(denovo_out.name, unit_output.name)

            with denovo_out.open() as denovo_f, unit_output.open() as unit_f:
                print(denovo_out.name)
                next(denovo_f), next(unit_f)
                old_scan = ""
                new_lines = dict()
                unit_lines = dict()
                for new_line in denovo_f:
                    new_scan = re.sub("(^[0-9]+\\:[0-9]+).*$", "\\1", new_line)
                    new_score = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\2", new_line)
                    if new_score == 'f\n':
                        new_score = 'inf'
                    new_score = float(new_score)
                    new_line = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\1\\3", new_line)

                    if old_scan == new_scan:
                        if new_line not in new_lines.keys():
                            new_lines[new_line] = new_score
                        elif new_lines[new_line] > new_score:
                            new_lines[new_line] = new_score
                    else:
                        unit_scan = old_scan
                        while (unit_scan == old_scan):
                            unit_line = unit_f.readline()
                            unit_scan = re.sub("(^[0-9]+\\:[0-9]+).*$", "\\1", unit_line)
                            unit_score = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\2", unit_line)
                            if unit_score == 'f\n':
                                unit_score = 'inf'
                            unit_score = float(unit_score)
                            unit_line = re.sub("(^[0-9]+\\:[0-9]+\t[^0-9]+)([-.inf0-9]+)(\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+$)", "\\1\\3", unit_line)

                            if unit_scan == old_scan:
                                if unit_line not in unit_lines.keys():
                                    unit_lines[unit_line] = unit_score
                                elif unit_lines[unit_line] > unit_score:
                                    unit_lines[unit_line] = unit_score
                        old_scan = new_scan

                        new_lines_conf = {k: new_lines[k] for k in new_lines if new_lines[k] < score_cutoff}
                        unit_lines_conf = {k: unit_lines[k] for k in unit_lines if unit_lines[k] < score_cutoff}

                        ## Checking that predictions with a good enough score are in both results.
                        self.assertTrue(set(new_lines_conf.keys()).issubset(unit_lines.keys()))
                        self.assertTrue(set(unit_lines_conf.keys()).issubset(new_lines.keys()))

                        ## Checking that score of good predictions is within no more than 0.05 when compared to unit test.
                        for k in unit_lines_conf.keys():
                            self.assertAlmostEqual((unit_lines_conf[k] - new_lines[k]) * 2, 0, places = 1)
                        for k in new_lines_conf.keys():
                            self.assertAlmostEqual((unit_lines[k] - new_lines_conf[k]) * 2, 0, places = 1)
                        new_lines = dict()
                        unit_lines = dict()
                        new_lines[new_line] = new_score
                        unit_lines[unit_line] = unit_score

        

if __name__ == '__main__':
    unittest.main()










