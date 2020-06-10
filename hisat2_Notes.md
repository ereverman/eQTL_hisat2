# Notes on my process for figuring out how to apply hisat to DSPR RNAseq data.

1. For troubleshooting, open an interactive node and create a conda environment:
```
srun --time=4:00:00 --ntasks=1 --nodes=1 --partition=sixhour --pty /bin/bash -l

module load anaconda
conda create -n HISAT2_env
# note that this may not make an environment that the final array script will be able to access. It will have a different name.
```
