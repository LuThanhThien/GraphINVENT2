import sys
from pathlib import Path

class DefaultConfig:
    """
    Configuration for running GraphINVENT2 jobs. Modify these parameters as necessary.
    """

    def update_paths(self, job_dir, tensorboard_dir):
        """
        Update dynamic paths for each job in the configuration.
        """
        self.params['job_dir'] = str(job_dir)+"/"
        self.params['tensorboard_dir'] = str(tensorboard_dir)+"/"

def load_config(args):
    """
    Load the configuration from a given path. This function is a placeholder and should be
    replaced with the actual implementation to load the configuration.

    Args:
        path (str): Path to the configuration file.

    Returns:
        config: Loaded configuration object.
    """
    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {cfg_path}")
    
    # Assert python file
    if not cfg_path.suffix == ".py":
        raise ValueError(f"Configuration file must be a Python file: {cfg_path}")
    
    # Load Config class from the specified path
    sys.path.append(str(cfg_path.parent))
    module_name = cfg_path.stem
    module = __import__(module_name)
    config_class = getattr(module, "Config")
    config = config_class()
    setattr(config, "python_path", args.python_path)
    return config

def dict2str(args: str):
    """
    Convert a dictionary to a string representation.

    Args:
        args (dict): Dictionary to convert.

    Returns:
        str: String representation of the dictionary.
    """
    return " ".join([f"--{k} {v}" for k, v in args.items()]) if isinstance(args, dict) else str(args)

def str2dict(args: str):
    """
    Convert a string representation of arguments back to a dictionary.

    Args:
        args (str): String representation of arguments.

    Returns:
        dict: Dictionary representation of the arguments.
    """
    if isinstance(args, dict):
        return args
    return {k: v for k, v in (item.split() for item in args.split("--") if item)}
