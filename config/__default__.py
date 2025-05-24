

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
