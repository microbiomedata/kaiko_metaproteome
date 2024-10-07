import os
import glob
import random
import subprocess
import yaml

from pathlib import Path, PureWindowsPath

def run_kaiko_denovo(config_path):
    kaiko_defaults_path = Path('kaiko_defaults.yaml')
    config = yaml.safe_load(kaiko_defaults_path.open())

    user_config_path = Path(config_path)
    assert user_config_path.exists(), "File " + str(user_config_path.absolute()) + " does not exist."
    config_user = yaml.safe_load(user_config_path.open())

    ## Overriding defaults if values found in user config.
    for section in config_user.keys():
        for param in config_user[section].keys():
            config[section][param] = config_user[section][param]

    mgf_dir = Path(PureWindowsPath(config['denovo']['mgf_dir']).as_posix())
    prefix = mgf_dir.name

    denovout_dir = Path('Kaiko_volume/Kaiko_intermediate/denovo_output/' + prefix)
    if not denovout_dir.exists():
        denovout_dir.mkdir()

    ## Step 1. Run Denovo using subprocess.
    kaiko_1_args = ["python", "src/kaiko_main.py", 
                    "--mgf_dir", mgf_dir.resolve(), 
                    "--train_dir", "model/",
                    "--decode_dir", denovout_dir.resolve(),
                    "--profile", config['denovo']['profile']]

    if config['denovo']['topk']:
        kaiko_1_args = kaiko_1_args + ["--topk"]
    if config['denovo']['multi_decode']:
        kaiko_1_args = kaiko_1_args + ["--multi_decode"]
    if config['denovo']['beam_search']:
        kaiko_1_args = kaiko_1_args + ["--beam_search", "--beam_size", config['denovo']['beam_size']]     

    print("DeNovo: Running the following command:\n")
    for i in range(len(kaiko_1_args)):
        kaiko_1_args[i] = str(kaiko_1_args[i])

    if not config['denovo']['cached']:
        print(" ".join(kaiko_1_args) + "\n")
        subprocess.run(kaiko_1_args, cwd = "Kaiko_denovo")
        # subprocess.call(kaiko_1_args, cwd = "Kaiko_denovo")

def run_new_parameters(topk, beam_size):
    config = {}

    prefix = os.path.split('Kaiko_volume/Kaiko_input_files/mgf_large_unit_test')[-1] 
    prefix = prefix + "_topk_{a1}_beam_size_{a4}".format(a1 = topk, a4 = beam_size)
    denovout_dir = 'Kaiko_volume/Kaiko_intermediate/denovo_output/' + prefix
    if not os.path.exists(denovout_dir):
        os.mkdir(denovout_dir)

    config['denovo'] = {'topk' : topk, 
                        'multi_decode' : True,
                        'beam_search' : True,
                        'beam_size' : beam_size,
                        'mgf_dir' : 'Kaiko_volume/Kaiko_input_files/mgf_large_unit_test',
                        'train_dir' : 'model/'}

    f = open(denovout_dir + '/config_unit_test.yaml', 'w')
    yaml.safe_dump(config, f)

    config_user = Path(denovout_dir + '/config_unit_test.yaml')
    denovout_dir = Path(denovout_dir)
    config_user = yaml.safe_load(config_user.open())

    for section in config_user.keys():
        for param in config_user[section].keys():
            config[section][param] = config_user[section][param]

    mgf_dir = Path(PureWindowsPath(config['denovo']['mgf_dir']).as_posix())
    prefix = mgf_dir.name

    ## Step 1. Run Denovo using subprocess.
    kaiko_1_args = ["python", "src/kaiko_main.py", 
                    "--mgf_dir", mgf_dir.resolve(), 
                    "--train_dir", "model/",
                    "--decode_dir", denovout_dir.resolve()]

    if config['denovo']['topk']:
        kaiko_1_args = kaiko_1_args + ["--topk"]
    if config['denovo']['multi_decode']:
        kaiko_1_args = kaiko_1_args + ["--multi_decode"]
    if config['denovo']['beam_search']:
        kaiko_1_args = kaiko_1_args + ["--beam_search", "--beam_size", config['denovo']['beam_size']]     

    print("DeNovo: Running the following command:\n")
    for i in range(len(kaiko_1_args)):
        kaiko_1_args[i] = str(kaiko_1_args[i])
    print(" ".join(kaiko_1_args) + "\n")

    subprocess.run(kaiko_1_args, cwd = "Kaiko_denovo")


def make_new_test_input(mgf_dir, proportion = 0.007, seed = 117):
    random.seed(seed)
    mgf_dir_name = os.path.split(mgf_dir)[-1] + "_unit_test"
    new_dir = os.path.join(os.path.dirname(mgf_dir), mgf_dir_name)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)

    mgf_files = glob.glob(mgf_dir + "/*.mgf")
    print(mgf_files)
    for mgf_file in mgf_files:
        new_lines = []
        with open(mgf_file) as file:
            for line in file:
                if line == "BEGIN IONS\n":
                    accept = random.uniform(0, 1)
                    print(accept)
                    print(len(new_lines))
                if accept < proportion:
                    new_lines += [line]

        mgf_test = os.path.split(os.path.splitext(mgf_file)[0])[-1]
        mgf_test = os.path.join(new_dir, "{mgf_test}_unit_test_prop_{proportion}_seed_{seed}.mgf".format(mgf_test = mgf_test, proportion = proportion, seed = seed))

        with open(mgf_test, 'a') as file:
            for line in new_lines:
                file.write(line)

    





