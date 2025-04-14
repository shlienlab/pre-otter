#### long-read demo-data usage ####
We have provided a demo dataset with expected output to help the user get started. 

The pipeline can be run with the demo data using `--test` as in the following command: 

```
singularity exec -e \ 
    -B /path/to/output/dir:/output \ 
    /singularity_cache/pre-otter.sif \
    pre-otter long-read \ 
    --test \
    -t my_API_token \ 
    -e email1@host.com,email2@host.com \ 
    -s  
```

Note, this is equivalent to running the following but on hg38 with fusion detection on chr20 only :
```
singularity exec -e \ 
    -B /path/to/output/dir:/output \ 
    /singularity_cache/pre-otter.sif \
    pre-otter long-read \ 
    -i /test_data/long_read/K562_ont_directRNA_10k_subsample.fastq.gz \ 
    -d \
    -p test \ 
    -f \
    -t my_API_token \ 
    -e email1@host.com,email2@host.com \ 
    -s  
```

The expected output from the test is also in the `/test_data` folder.

Note, the demo dataset is a small subset of data from [Chen et al. Nature Methods (2025)](https://www.nature.com/articles/s41592-025-02623-4). For more information, check the [SG-NEx](https://github.com/GoekeLab/sg-nex-data/) website. 