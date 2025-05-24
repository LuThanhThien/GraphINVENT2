"""
Example submission script for a GraphINVENT2 training job (unconditional generation,
 not fine-tuning/optimization job). This can be used to pre-train a model before
 a reinforcement learning (fine-tuning) job.

To run, type:
  user@cluster GraphINVENT2$ python submit.py config_path [--python_path python_path]
"""
# load general packages and functions
import csv
import sys
import os
from pathlib import Path
import subprocess
import time
from config import load_config, parse_args_submit

def submit(config):
    """
    Creates and submits submission scripts based on the provided configuration.

    Args:
        config: Configuration object containing all settings and parameters.
    """
    if config.job_type == "pre-filter" or config.job_type == "post-filter":
        submit_filter(config)
        return
    
    dataset_output_path, tensorboard_path = create_output_directories(config)
    submit_jobs(config, dataset_output_path, tensorboard_path)

def submit_one_filter(config, params, script_name):
    """
    Creates and submits a single job script based on the provided parameters.

    Args:
        config: Configuration object containing all settings and parameters.
        params: Dictionary of parameters for the job.
        script_name: Name of the script to be executed.
    """
    # This path as PYTHONPATH
    proj_path = os.path.dirname(os.path.abspath(__file__))
    # Export the project path to PYTHONPATH
    os.environ["PYTHONPATH"] = proj_path

    process_run = [
        config.python_path, 
        script_name,
    ]
    for key, value in params.items():
        process_run.append(f"--{key}")
        process_run.append(str(value))

    subprocess.run(process_run, check=True)

def submit_filter(config):
    """
    Creates and submits submission scripts based on the provided configuration.

    Args:
        config: Configuration object containing all settings and parameters.
    """
    submit_one_filter(config, config.filter_params, "filter/filter_smi.py")
    
def create_output_directories(config):
    """
    Creates output and tensorboard directories.

    Args:
        config: Configuration object containing dataset and jobname.

    Returns:
        A tuple of dataset_output_path and tensorboard_path.
    """
    base_path = Path(f"./output/output_{config.dataset}")
    # dataset_output_path = base_path / config.jobname if config.jobname else base_path
    dataset_output_path = base_path 
    tensorboard_path = dataset_output_path / "tensorboard"

    dataset_output_path.mkdir(parents=True, exist_ok=True)
    tensorboard_path.mkdir(parents=True, exist_ok=True)
    print(f"* Creating dataset directory {dataset_output_path}/", flush=True)

    return dataset_output_path, tensorboard_path

def submit_jobs(config, dataset_output_path, tensorboard_path):
    """
    Submits the specified number of jobs by creating subdirectories and submission scripts.

    Args:
        config: Configuration object containing job and execution details.
        dataset_output_path: Path object for the dataset output directory.
        tensorboard_path: Path object for the tensorboard directory.
    """
    jobdir_end_idx = config.jobdir_start_idx + config.n_jobs
    for job_idx in range(config.jobdir_start_idx, jobdir_end_idx):
        job_dir = dataset_output_path / f"job_{job_idx}/"
        tensorboard_dir = tensorboard_path / f"job_{job_idx}/"

        tensorboard_dir.mkdir(parents=True, exist_ok=True)
        create_job_directory(job_dir, config)

        config.update_paths(job_dir, tensorboard_dir)
        write_input_csv(config, job_dir, filename="input.csv")
        submit_job(config, job_dir)

        print("-- Sleeping 2 seconds.")
        time.sleep(2)

def create_job_directory(job_dir, config):
    """
    Creates a job directory, handles overwriting based on configuration.

    Args:
        job_dir: Path object for the job directory.
        config: Configuration object with job type and overwrite settings.
    """
    try:
        job_dir.mkdir(parents=True, exist_ok=config.force_overwrite or config.job_type in ["generate", "test"])
        print(f"* Creating model subdirectory {job_dir}/", flush=True)
    except FileExistsError:
        print(f"-- Model subdirectory {job_dir} already exists.", flush=True)
        if not config.restart:
            return

def submit_job(config, job_dir):
    """
    Writes and submits a job based on the SLURM configuration or runs directly.

    Args:
        config: Configuration object with paths and SLURM settings.
        job_dir: Path object for the job directory.
    """
    if config.use_slurm:
        print("* Writing submission script.", flush=True)
        write_submission_script(config, job_dir)

        print("* Submitting job to SLURM.", flush=True)
        subprocess.run(["sbatch", str(job_dir / "submit.sh")], check=True)
    else:
        print("* Running job as a normal process.", flush=True)
        # Get python path from current environment
        subprocess.run([config.python_path, config.graphinvent_path + "main.py", "--job-dir", str(job_dir)+"/"], check=True)

def write_input_csv(config, job_dir, filename="params.csv") -> None:
    """
    Writes job parameters/hyperparameters from the config object to a CSV file.
    Args:
        config: Configuration object containing all settings and parameters.
        job_dir (Path): The directory where the job will run.
        filename (str): Filename for the CSV output, default is "params.csv".
    """
    dict_path = job_dir / filename

    try:
        # open the file at dict_path in write mode
        with dict_path.open(mode="w", newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter=";")
            for key, value in config.params.items():
                writer.writerow([key, value])
    except IOError as e:
        # exception handling for any I/O errors
        print(f"Failed to write file {dict_path}: {e}")

def write_submission_script(config, job_dir) -> None:
    """
    Writes a submission script (`submit.sh`) using settings from the config object.
    Args:
        config: Configuration object containing all settings and parameters.
        job_dir (Path): The directory where the job will run.
    """
    submit_filename = job_dir / "submit.sh"
    output_filename = job_dir / f"output.o${{SLURM_JOB_ID}}"
    main_py_path = Path(config.graphinvent_path) / "main.py"

    with submit_filename.open("w") as submit_file:
        submit_file.write("#!/bin/bash\n")
        submit_file.write(f"#SBATCH -A {config.account}\n")
        submit_file.write(f"#SBATCH --job-name={config.job_type}\n")
        submit_file.write(f"#SBATCH --time={config.run_time}\n")
        submit_file.write("#SBATCH --gpus-per-node=T4:1\n")
        submit_file.write("hostname\n")
        submit_file.write("export QT_QPA_PLATFORM='offscreen'\n")
        submit_file.write(f"({config.python_path} {main_py_path} --job-dir {job_dir} > {output_filename})\n")

if __name__ == "__main__":
    args = parse_args_submit()
    config = load_config(args)  # load the config object
    submit(config)     # pass the config object to the submit function
