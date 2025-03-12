#### short-read demo-data usage ####
We have provided a demo dataset with expected output to help the user get started. The pipeline can be run with the demo data using the following command: 
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif \
    pre-otter short-read \
    -r1 /test_data/short_read/K562_Illumina_replicate4_200K_subsample_R1.fastq.gz \
    -r2 /test_data/short_read/K562_Illumina_replicate4_200K_subsample_R2.fastq.gz \
    -p pre-otter_short_read_test \
    -t $MY_OTTER_API_TOKEN \
    -e email1@host.com,email2@host.com \
    -S  
```
The expected output from the test is also in the `/test_data` folder.