from config.__default import DefaultConfig

class Config(DefaultConfig):
    """
    Configuration for running GraphINVENT2 jobs. Modify these parameters as necessary.
    """
    def __init__(self):
        # set paths here
        self.graphinvent_path = "./graphinvent/"
        self.data_path        = "./data/pre-training/"
        
        # set SLURM params here (if using SLURM)
        self.use_slurm        = False              # use SLURM or not
        self.run_time         = "0-06:00:00"       # d-hh:mm:ss
        self.account          = "XXXXXXXXXX"       # if cluster requires specific allocation/account, use here

        # define what you want to do for the specified job(s)
        self.dataset          = "minicoconut"  # dataset name in "./data/pre-training/"
        self.job_type         = "pre-filter"   # "pre-filter", " post-filter", "preprocess", "train", "generate", "fine-tune", or "test"
        self.jobdir_start_idx = 0              # where to start indexing job dirs
        self.n_jobs           = 1              # number of jobs to run per model
        self.restart          = False          # whether or not this is a restart job
        self.force_overwrite  = True           # overwrite job directories which already exist
        self.jobname          = self.job_type  # label used to create a job sub directory (can be anything)

        # define dataset-specific parameters
        self.params = {
            "atom_types"     : ["C", "N", "O", "S"],              
            "formal_charge"  : [-1, 0, +1,],                                         
            "max_n_nodes"    : 115,                                        
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

        # config for pre-filter
        self.filter_params = {
            "input_file": f"./data/pre-training/minicoconut/raw/flavonoid_new_filter.smi",
            "output_file": f"./data/pre-training/minicoconut/raw/flavonoid_final.smi",
            "overwrite": True,
            "derivative": {
                "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c21",  # naringenin
                "O=c1ccc2ccccc2o1",      # flavone
                "O=c1cc(-c2ccc(O)cc2)oc2ccc(O)cc12",   # apigenin
                "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)c(O)cc12",  # quercetin
                "O=c1cc(-c2ccc(O)c(O)c2)oc2ccc(O)cc12",     # kaempferol
                # "O=c1cc(-c2ccc(O)c(O)c2)oc2ccc(O)cc12",    # luteolin
                "O=c1cc(-c2ccccc2)oc2ccccc12",              # chrysin
                # "O=c1cc(-c2ccc(O)c(O)c2)oc2ccc(O)cc12",     # galangin
                # "O=c1cc(-c2ccc(O)c(O)c2)oc2ccc(O)cc12",     # fisetin
                # "O=c1cc(-c2ccc(O)c(O)c2)oc2ccc(O)cc12",      # baicalein
                "O=C1c2c(O)cc(O)cc2OC(c2ccc(O)c(O)c2)C1O",   #taxifolin
                "O=c1c(O)c(-c2ccc(O)cc2O)oc2cc(O)cc(O)c12"   #Morin
            },
            # "included_atoms": {'C', 'N', 'O', 'S'},
            # "included_charges": {-1, 0, 1},
            # "max_n_nodes": 50,
            "min_n_nodes": 10,
            # "max_n_molecules": 20000,
        }

        # config for post filter
        #TODO: Add postprocess filter parameters if needed
