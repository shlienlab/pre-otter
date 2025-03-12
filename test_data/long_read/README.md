#### long-read demo-data usage ####
We have provided a demo dataset with expected output to help the user get started. The pipeline can be run with the demo data using the following command: 
```
singularity exec -e \ 
    -B /path/to/output/dir:/output \ 
    /singularity_cache/pre-otter.sif \
    pre-otter long-read \ 
    -i /test_data/long_read/K562_ont_directRNA_10k_subsample.fastq.gz \ 
    -d 
    -p test \ 
    -t my_API_token \ 
    -e email1@host.com,email2@host.com \ 
    -s  
```
The expected output from the test is also in the `/test_data` folder.