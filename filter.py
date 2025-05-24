"""
Example submission script for a GraphINVENT2 training job (unconditional generation,
 not fine-tuning/optimization job). This can be used to pre-train a model before
 a reinforcement learning (fine-tuning) job.

To run, type:
 user@cluster GraphINVENT2$ python submit.py config_path [--python_path python_path]
"""
# load general packages and functions
import sys
import os
from argparse import ArgumentParser
from config import load_config
import subprocess

def parse_args():
    """
    Parse command line arguments for the script.
    """
    default_python_path = os.path.join(os.path.dirname(sys.executable), "python")
    parser = ArgumentParser(description="Submit GraphINVENT2 filter jobs.")
    parser.add_argument("config", type=str, help="Path to the configuration file.")
    parser.add_argument("--python_path", type=str, default=default_python_path, help="Path to the Python executable.")
    args = parser.parse_args()
    return args

def submit(config):
    """
    Creates and submits submission scripts based on the provided configuration.

    Args:
        config: Configuration object containing all settings and parameters.
    """
    subprocess.run([config.python_path, config.tool_path + "filter/filter_smi.py"], check=True)

if __name__ == "__main__":
    args = parse_args()
    config = load_config(args)  # load the config object
    submit(config)     # pass the config object to the submit function
