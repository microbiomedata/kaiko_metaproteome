import argparse

from unit_test_util import run_kaiko_denovo


print("Running just denovo inference on the given input.\n")
parser = argparse.ArgumentParser()
parser.add_argument("-config", "--config", help="YAML file with parameters for pipeline, and input files path. Any found values in config replace the default values.")


args = parser.parse_args()
assert args.config is not None, "Please provide a config file with the spectrum path (mgf files). Use flag --config to pass the path. See the tutorial for more information."

run_kaiko_denovo(args.config)
