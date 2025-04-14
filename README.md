# Welcome to pre-otter #

### Rational ###

* This repository is for providing a tool for uniform and reproducible gene expression interoperable with the otter web app.
* The pre-built container can be found here: https://hub.docker.com/repository/docker/shlienteam/pre-otter

### How do I get set up? ###
1. Have at least 48GB of memory, but certain containerization software could require more to run.
2. The more CPUs the merrier, we tested with 16 and 24 CPUs.
3. Singularity must be installed with `/tmp` being applied, or a directory must be bound to the `/tmp` container directory.
4. Have enough free space in `/tmp` for the references, index, and bam files to be created transiently.

#### Singularity ####
```
singularity pull /singularity_cache/pre-otter.sif docker://shlienteam/pre-otter:1.1.0
```
Symlinked paths are not supported. Singularity currently cannot automatically resolve home directories located on symlinked paths. As a result, the home directory must be explicitly bound using a bind path, rather than relying on the automatic mounting of the home directory.

### Usages: ###
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter -h
    
usage: pre-otter [-h] {short-read,long-read} ...

Process RNA-seq fastqs to gene expression counts, TPM.

positional arguments:
  {short-read,long-read}
                        Mode of operation
    short-read          Run in short-read mode.
    long-read           Run in long-read mode.

options:
  -h, --help            show this help message and exit
```

#### short-read usage: ####
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter short-read -h

usage: pre-otter short-read [-h] [-r1] [-r2] [-p] [--om] [--v2] [-t] [-e] [-s]

options:
  -h  , --help     show this help message and exit
  -r1 , --read1    fastq read1 files, comma separated if more than one.
  -r2 , --read2    fastq read2 files, comma separated if more than one.
  -p  , --prefix   prefix to be used for the output files.
  --om             use otter model instead of hierarchical model.
  --v2             use version 2 of model instead of version 1.
  -t  , --token    otter web app API token.
  -e  , --email    comma separate list of email addresses, otter web app will email when analysis complete.
  -s  , --save     NOTE: setting this will allow the otter web app to keep the TPM file sent.
```

Using a single pair of fastq files:
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter short-read \
    -r1 /data/read1_R1.fastq.gz \
    -r2 /data/read2_R2.fastq.gz \
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
    /singularity_cache/pre-otter.sif pre-otter short-read \
    -r1 /data/xyz_Lane1_R1.fastq.gz,/data/xyz_Lane2_R1.fastq.gz,/data/xyz_Lane3_R1.fastq.gz \
    -r2 /data/xyz_Lane1_R2.fastq.gz,/data/xyz_Lane2_R2.fastq.gz,/data/xyz_Lane3_R2.fastq.gz \
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

cd /path/to/run/directory/pre-otter-run

singularity exec -e \
-B /directory/path/to/store/reference/indexes:/reference \
-B /path/to/fastqs/directory/:/data \
-B /path/where/you/want/output/stored/:/output \
/singularity_cache/pre-otter.sif pre-otter short-read \
-r1 /data/SRR30668570_1.fastq.gz \
-r2 /data/SRR30668570_2.fastq.gz \
-p SRR30668570 \
-t e01918ac92401___Your_Otter_Web_App_API_TOKEN___1ef8800fa163e7415e4 \
-s
```

#### long-read usage: ####
```
usage: pre-otter long-read [-h] [-i] [-c] [-d] [-f] [--test] [-p] [--om] [--v2] [-t] [-e] [-s]

options:
  -h , --help     show this help message and exit
  -i , --input    input file(s) seperated by comma or single input directory.
  -c , --cdna     cDNA-PCR library kit used, use PCB114 by default.
  -d , --drna     direct-RNA library specified, skipping read trimming, set to FALSE by default.
  -f, --fusion    fusion gene detection, set to FALSE by default.
  --test          run a test using demo data with predefined options.
  -p , --prefix   prefix to be used for the output files.
  --om            use otter model instead of hierarchical model.
  --v2            use version 2 of model instead of version 1.
  -t , --token    otter web app API token.
  -e , --email    comma separate list of email addresses, otter web app will email when analysis complete.
  -s , --save     NOTE: setting this will allow the otter web app to keep the TPM file sent.
```

Using cDNA-PCR library:
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter long-read \
    -i ont_read.fastq.gz \
    -c PCB114 \
    -f \
    -p my_counts_prefix \
    -t my_API_token \
    -e email1@host.com,email2@host.com \
    -s
```

Using direct-RNA library:
```
singularity exec -e \
    -B /path/to/reference/dir:/reference \
    -B /path/to/data/dir:/data \
    -B /path/to/output/dir:/output \
    /singularity_cache/pre-otter.sif pre-otter long-read \
    -i ont_read.fastq.gz \
    -d \
    -f \
    -p my_counts_prefix \
    -t my_API_token \
    -e email1@host.com,email2@host.com \
    -s
```

#### demo-data usage: ####
For instructions to test the provided demo datasets, please see the README found under [test_data/short_read](test_data/short_read/README.md) and [test_data/long_read](test_data/long_read/README.md).

#### Additional Notes ####
- Long-read usage supports cDNA kits in any form of `{PCS109,PCS110,PCS111,PCS114,LSK114,PCB111,PCB114}`. 
- Long-read usage accepts directory input including the entire binded folder: `-i /data` or a subdirectory folder: `-i /data/subdir`.
- For fusion detection (`-f`), the genome index reference (~14 GB) will be downloaded to the binded `/reference` folder, if not detected. To save time, the user can download and extract the reference files from [here](https://figshare.com/ndownloader/files/25410494) with the download [script](scripts/download_jaffal_references.sh).
- In the case of no gene fusion detected, pre-otter will output an empty fusion results file.
- Please allow up to 3-4 hours for the build process to complete. We appreciate your patience during this time.

### Who do I talk to? ###
```
Contact:
- Scott Davidson <scott.davidson@sickkids.ca>
- Pedro Lemos Ballester <pedro.lemosballester@sickkids.ca>
- Sandy Fong <sandy.fong@sickkids.ca>
```