# Welcome to pre-otter #

### Rational ###

* This repository is for providing a tool for uniform and reproducible gene expression interoperable with the otter web app.

### How do I get set up? ###
1. Have at least 48GB to 60GB of memory.
2. The more cores the merrier, we tested with 16 cores.
3. Singularity must be installed with `/tmp` being applied by admin, or a directory must be bound to the `/tmp` container directory.
4. Have enough free space in `/tmp` for the references, index, and bam files to be created transiently.

#### Singularity ####
```
singularity pull /singularity_cache/pre-otter.sif docker://shlienteam/pre-otter:0.0.1
```

Usages:
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter -h
    
usage: pre-otter [-h] [-r1] [-r2] [-p] [-t] [-e] [-s]

Process RNA-seq fastqs to gene expression counts, TPM.

options:
  -h, --help      show this help message and exit
  -r1 , --read1   fastq read1 files, comma separated if more than one pair.
  -r2 , --read2   fastq read2 files, comma separated if more than one pair.
  -p , --prefix   output prefix.
  -t , --token    Otter web app API token.
  -e , --email    Comma separate list of email addresses, otter web app will
                  email when analysis complete.
  -s, --save      NOTE: Setting this will allow the otter web app to keep the
                  TPM file sent.
```
Using a single pair of fastq files:
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter \
    -r1 read1_R1.fastq.gz \
    -r2 read2_R2.fastq.gz \
    -p my_counts_prefix \
    -t my_API_token \
    -e email1@host.com,email2@host.com \
    -s
```

Using fastq files for the same sample sequenced across multiple lanes:
```
singularity exec -e 
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter \
    -r1 xyz_Lane1_R1.fastq.gz,xyz_Lane2_R1.fastq.gz,xyz_Lane3_R1.fastq.gz \
    -r2 xyz_Lane1_R2.fastq.gz,xyz_Lane2_R2.fastq.gz,xyz_Lane3_R2.fastq.gz \
    -p /path/to/my/output/directory/my_counts_prefix
    -t my_API_token \
    -e email1@host.com \
    -s
```
#### On HPC with Slurm Scheduler ####
```
#!/bin/bash
#SBATCH --job-name pre_otter
#SBATCH -e /path/to/logs/directory/%x.e%j.out
#SBATCH -o /path/to/logs/directory/%x.o%j.out
#SBATCH --get-user-env
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --mem 48G
#SBATCH -t 24:00:00
#SBATCH -N 1 -c 8


module load Singularity/3.11.3

cd /hpf/largeprojects/adam/projects/pedro/pre-otter-run

singularity exec -e \
-B /directory/path/to/store/reference/indexes:/reference \
-B /path/to/fastqs/directory/:/data \
-B /path/where/you/want/output/stored/:/output \
/singularity_cache/pre-otter.sif pre-otter \
-r1 /data/SRR30668570_1.fastq.gz \
-r2 /data/SRR30668570_2.fastq.gz \
-p SRR30668570 \
-t e01918ac92401___Your_Otter_Web_App_API_TOKEN___1ef8800fa163e7415e4
-s
```

### Who do I talk to? ###
```
Contact:
- Scott Davidson <scott.davidson@sickkids.ca>
- Pedro Lemos Ballester <pedro.lemosballester@sickkids.ca>
```