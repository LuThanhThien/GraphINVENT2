from __default__ import DefaultConfig

class Config(DefaultConfig):
    """
    Configuration for running GraphINVENT2 jobs. Modify these parameters as necessary.
    """
    def __init__(self):
        # set paths here
        self.graphinvent_path = "./graphinvent/"
        self.data_path        = "./data/pre-training/"
        
        # set SLURM params here (if using SLURM)
        self.use_slurm        = False               # use SLURM or not
        self.run_time         = "0-06:00:00"       # d-hh:mm:ss
        self.account          = "XXXXXXXXXX"       # if cluster requires specific allocation/account, use here

        # define what you want to do for the specified job(s)
        self.dataset          = "coconut"      # dataset name in "./data/pre-training/"
        self.job_type         = "preprocess"   # "preprocess", "train", "generate", "fine-tune", or "test"         
        self.jobdir_start_idx = 0              # where to start indexing job dirs
        self.n_jobs           = 1              # number of jobs to run per model
        self.restart          = False          # whether or not this is a restart job
        self.force_overwrite  = True           # overwrite job directories which already exist
        self.jobname          = self.job_type  # label used to create a job sub directory (can be anything)
        
        # define dataset-specific parameters
        self.params = {
            "atom_types"     : ["C", "N", "O", "S"],              
            "formal_charge"  : [-1, 0, +1,],                                         
            "max_n_nodes"    : 120,                                        
            "job_type"       : self.job_type,
            "dataset_dir"    : f"{self.data_path}{self.dataset}/",
            "restart"        : self.restart,
            "sample_every"   : 1,
            "init_lr"        : 1e-4,
            "epochs"         : 100,
            "batch_size"     : 500,
            "block_size"     : 10000,
            "device"         : "cuda",  # or "cpu" if no CUDA
            "n_samples"      : 1000, 
            "generation_epoch": 100,  #** <-- which model to use (i.e. which epoch)
            # additional paramaters can be defined here, if different from the "defaults"
            # for instance, for "generate" jobs, don't forget to specify "generation_epoch"
            # and "n_samples"
        }
