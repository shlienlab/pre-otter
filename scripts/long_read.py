#!/usr/bin/env python3

import logging
import os
import subprocess
import argparse
import configparser
import tempfile
import shutil
import multiprocessing
import json
import pandas as pd
import requests as reqs
import time

from typing import Any
    

def load_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    return config['defaults']  # Return the default values from the config file


def validate_input_files(input_path: str):
    """
    Check if the input is a directory or a comma-separated list of files.
    Validate that all files are either .fastq or .fastq.gz.
    :param input_path: (str) Path inputted by user
    """
    valid_files = []  # List to hold valid files
    total_count = 0  # Counter for total valid .fastq files

    if len(input_path) < 2 and os.path.isdir(input_path[0]):
        # Check the root directory itself for valid files
        root_files = os.listdir(input_path[0])
        valid_root_files = [f for f in root_files if f.endswith(('.fastq', '.fastq.gz'))]
        
        # If any valid files are found in the root directory, add to the list and count
        if valid_root_files:
            logging.error(f"Valid root directory found: {input_path[0]}")
            valid_files.extend([os.path.join(input_path[0], f) for f in valid_root_files])
            total_count += len(valid_root_files)

        # Recursively search through subdirectories
        for root, dirs, files in os.walk(input_path[0]):
            for dir_name in dirs:
                sub_dir_path = os.path.join(root, dir_name)
                sub_files = os.listdir(sub_dir_path)
                
                # Check if any valid .fastq or .fastq.gz files exist in the subdirectory
                valid_sub_files = [f for f in sub_files if f.endswith(('.fastq', '.fastq.gz'))]
                if valid_sub_files:  # If there are valid files in the subdirectory
                    logging.info(f"Valid subdirectory found: {sub_dir_path}")
                    valid_files.extend([os.path.join(sub_dir_path, f) for f in valid_sub_files])
                    total_count += len(valid_sub_files)

        # Return the list of file paths
        logging.info(f"Total count of read files inputted: {total_count}")
        if valid_files:
            return valid_files
        else:
            raise ValueError(f"No valid .fastq or .fastq.gz files found in {input_path}")

    # Check if input is a comma-separated list of files
    else:
        for file in input_path:
            file = file.strip()  # Remove any surrounding whitespace
            # Check if the file ends with .fastq or .fastq.gz
            if not (file.endswith('.fastq') or file.endswith('.fastq.gz')):
                raise ValueError(f"Invalid input: {file}. Input must be a valid directory or file with extension: .fastq or .fastq.gz.")
            
            # Check if "/data" + file_name exists
            if file.startswith("/data/") or file.startswith("/test_data/"):
                data_path = file
            else:
                data_path = os.path.join("/data/", file)  # Construct the full path to /data/file_name
                if not os.path.exists(data_path):
                    data_path = os.path.join("/test_data/", file)
            
            if not os.path.exists(data_path):
                raise ValueError(f"File not found: {data_path}. Please ensure it exists.")
            valid_files.append(data_path)
        
        logging.info(f"Total count of read files inputted: {len(valid_files)}")
        return valid_files

    
def concatenation(tmp_dir: str, prefix: str, input: list, fastcat_opts: str):
    """
    Performs concatenation of multiple input FASTQ files using FastCat
    :param tmp_dir: (str) Temporary directory for storing intermediate files
    :param input_dir: (str) List of .FASTQ files inputted by user
    :param fastcat_opts: (str) Options for fastcat
    """
    # fastcat require an uncreated folder for --histogram
    # Define the output directories and filenames
    fastcat_output_dir = f"{tmp_dir}/fastcat"
    os.makedirs(fastcat_output_dir, exist_ok=True)
    
    fastcat_output = f"{fastcat_output_dir}/{prefix}_seqs.fastq.gz"
    histograms_output = f"{fastcat_output_dir}/{prefix}_histograms"
    
    cpu_count = multiprocessing.cpu_count()
    # Run the fastcat command and pipe it to bgzip
    fastcat_command = [
        "fastcat", str(fastcat_opts), "-s", prefix, "--histograms", "histograms", *input
    ]
    logging.info(fastcat_command)
    with subprocess.Popen(fastcat_command, stdout=subprocess.PIPE) as proc:
        # Pipe the output of fastcat to bgzip
        with subprocess.Popen(
            ["bgzip", "-@", str(cpu_count)], stdin=proc.stdout, stdout=open(fastcat_output, 'wb')
        ) as bgzip_proc:
            bgzip_proc.communicate()  # Wait for bgzip to finish

    # Move histograms file
    subprocess.run(["mv", "histograms", histograms_output], check=True)

    logging.info(f"Output saved to: {fastcat_output}")
    logging.info(f"Histograms moved to: {histograms_output}")


def adaptor_trimming(tmp_dir: str, used_fastcat: bool, input: str, cdna_kit: str, pychopper_backend: str, prefix: str):
    """
    Performs adaptor trimming and full-length read extraction using PyChopper
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :param used_fastcat: (bool) Flag to use output generated from fastcat
    :param input: (str) If only 1 file is inputted, can directly use for adaptor triming
    :param cdna_kit: (str) cDNA kit {PCS109,PCS110,PCS111,PCS114,LSK114,PCB111,PCB114}
    :param pychopper_backend: (str) pychopper option for detection method.
    :param prefix: (str) The prefix to use for output files.
    """
    # Define the output directories and filenames
    if used_fastcat is True:
        file_input = f"{tmp_dir}/fastcat/{prefix}_seqs.fastq.gz"
    else:
        file_input = input
    
    cpu_count = multiprocessing.cpu_count()
    pychopper_output_dir = f"{tmp_dir}/pychopper"
    os.makedirs(pychopper_output_dir, exist_ok=True)
    
    pychopper_full_length_output = f"{pychopper_output_dir}/{prefix}_full_length_reads.fastq"
    
    # Run pychopper command
    pychopper_command = [
        "pychopper",
        "-t", str(cpu_count),
        "-k", cdna_kit,
        "-m", pychopper_backend,
        file_input,
        pychopper_full_length_output
    ]
    logging.info(pychopper_command)
    subprocess.run(pychopper_command, check=True)

    # Run bgzip command to compress the full-length reads
    bgzip_command = ["bgzip", "-@", str(cpu_count), pychopper_full_length_output]
    subprocess.run(bgzip_command, check=True)

    # Move pychopper.* files to the appropriate directory
    # subprocess.run(["mv", "pychopper.pdf", pychopper_output_dir], check=True)
    subprocess.run(f"mv pychopper.* {pychopper_output_dir}", shell=True, check=True)

    logging.info(f"Full-length reads saved to: {pychopper_full_length_output}.gz")
    logging.info(f"pychopper.* files moved to: {pychopper_output_dir}")


def reference_indexing(tmp_dir: str, minimap2_index_opts: str, ref_transcriptome: str):
    """
    Performs reference indexing using minimap2
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :param minimap2_index_opts: (str) Options for minimap2
    :param ref_transcriptome: (str) The path to the reference transcriptome file.
    """
    minimap2_output_dir = f"{tmp_dir}/minimap2"
    cpu_count = multiprocessing.cpu_count()
    os.makedirs(minimap2_output_dir, exist_ok=True)
    reference_index_output = f"{minimap2_output_dir}/ref.mmi"

    minimap2_command = [
        "minimap2",
        "-t", str(cpu_count),
        *minimap2_index_opts.replace('"', '').split(),
        "-d", reference_index_output,
        ref_transcriptome
    ]
    logging.info(minimap2_command)
    subprocess.run(minimap2_command, check=True)
    logging.info(f"Reference transcriptome index saved to: {reference_index_output}")


def drna_transcriptome_mapping(tmp_dir: str, used_fastcat: bool, input: str, minimap2_drna_tx_opts: str, aln_filtering_opts: str, prefix: str):
    """
    Performs dRNA transcriptome mapping using minimap2
    and creates sorted BAM file using samtools
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :param used_fastcat: (bool) Flag to use output generated from fastcat
    :param input: (str) If only 1 file is inputted, can directly use for mapping
    :param minimap2_cdna_tx_opts: (str) Options for minimap2
    :param ref_transcriptome: (str) The path to the reference transcriptome file.
    :param prefix: (str) The prefix to use for output files.
    """
    # Define the output directories and filenames
    if used_fastcat is True:
        file_input = f"{tmp_dir}/fastcat/{prefix}_seqs.fastq.gz"
    else:
        file_input = input

    minimap2_output_dir = f"{tmp_dir}/minimap2"
    cpu_count = multiprocessing.cpu_count()
    os.makedirs(minimap2_output_dir, exist_ok=True)

    # Checking reference transcriptome index file
    reference_index = f"{minimap2_output_dir}/ref.mmi"
    if not os.path.exists(reference_index):
        raise FileNotFoundError(f"Reference transcriptome file {reference_index} does not exist.")

    # Step 1: Run minimap2 and pipe it to samtools view
    aln_bam_output = f"{minimap2_output_dir}/{prefix}_aln.bam"
    minimap2_command = [
        "minimap2",
        "-t", str(cpu_count),
        *minimap2_drna_tx_opts.replace('"', '').split(),
        reference_index,
        file_input
    ]
    logging.info(minimap2_command)
    with subprocess.Popen(minimap2_command, stdout=subprocess.PIPE) as minimap2_proc:
        samtools_view_command = [
            "samtools", "view", "-Sb", *aln_filtering_opts.replace('"', '').split()
        ]
        logging.info(samtools_view_command)
        with subprocess.Popen(samtools_view_command, stdin=minimap2_proc.stdout, stdout=open(aln_bam_output, 'wb')) as samtools_proc:
            samtools_proc.communicate()  # Wait for samtools to finish

    # Step 2: Run samtools sort to sort the alignment BAM file
    sorted_bam_output = f"{minimap2_output_dir}/{prefix}_aln_sorted.bam"
    samtools_sort_command = [
        "samtools", "sort", "-@", str(cpu_count), aln_bam_output, "-o", sorted_bam_output
    ]
    logging.info(samtools_sort_command)
    subprocess.run(samtools_sort_command, check=True)

    logging.info(f"Aligned BAM file saved to: {aln_bam_output}")
    logging.info(f"Sorted BAM file saved to: {sorted_bam_output}")


def cdna_transcriptome_mapping(tmp_dir: str, minimap2_cdna_tx_opts: str, aln_filtering_opts: str, prefix: str):
    """
    Performs cDNA transcriptome mapping using minimap2
    and creates sorted BAM file using samtools
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :param minimap2_cdna_tx_opts: (str) Options for minimap2
    :param ref_transcriptome: (str) The path to the reference transcriptome file.
    :param prefix: (str) The prefix to use for output files.
    """
    # Define the output directories and filenames
    full_length_reads = f"{tmp_dir}/pychopper/{prefix}_full_length_reads.fastq.gz"
    minimap2_output_dir = f"{tmp_dir}/minimap2"
    cpu_count = multiprocessing.cpu_count()
    os.makedirs(minimap2_output_dir, exist_ok=True)

    # Checking reference transcriptome index file
    reference_index = f"{minimap2_output_dir}/ref.mmi"
    if not os.path.exists(reference_index):
        raise FileNotFoundError(f"Reference transcriptome file {reference_index} does not exist.")
    
    # Step 1: Run minimap2 and pipe it to samtools view
    aln_bam_output = f"{minimap2_output_dir}/{prefix}_aln.bam"
    minimap2_command = [
        "minimap2",
        "-t", str(cpu_count),
        *minimap2_cdna_tx_opts.replace('"', '').split(),
        reference_index,
        full_length_reads
    ]
    logging.info(minimap2_command)
    with subprocess.Popen(minimap2_command, stdout=subprocess.PIPE) as minimap2_proc:
        samtools_view_command = [
            "samtools", "view", "-Sb", *aln_filtering_opts.replace('"', '').split()
        ]
        logging.info(samtools_view_command)
        with subprocess.Popen(samtools_view_command, stdin=minimap2_proc.stdout, stdout=open(aln_bam_output, 'wb')) as samtools_proc:
            samtools_proc.communicate()  # Wait for samtools to finish

    # Step 2: Run samtools sort to sort the alignment BAM file
    sorted_bam_output = f"{minimap2_output_dir}/{prefix}_aln_sorted.bam"
    samtools_sort_command = [
        "samtools", "sort", "-@", str(cpu_count), aln_bam_output, "-o", sorted_bam_output
    ]
    logging.info(samtools_sort_command)
    subprocess.run(samtools_sort_command, check=True)

    logging.info(f"Aligned BAM file saved to: {aln_bam_output}")
    logging.info(f"Sorted BAM file saved to: {sorted_bam_output}")


def quantitation(tmp_dir: str, ref_transcriptome: str, prefix: str):
    """
    Performs transcript abundance quantification using Salmon
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :param ref_transcriptome: (str) The path to the reference transcriptome file.
    :param prefix: (str) The prefix to use for output files.
    """
    # Define output directories and filenames
    salmon_output_dir = f"{tmp_dir}/salmon/{prefix}"
    cpu_count = multiprocessing.cpu_count()
    os.makedirs(salmon_output_dir, exist_ok=True)

    # Step 1: Run salmon quant
    salmon_command = [
        "salmon", "quant",
        "-p", str(cpu_count),
        "--noErrorModel", 
        "-t", ref_transcriptome,
        "-l", "SF",
        "-a", f"{tmp_dir}/minimap2/{prefix}_aln_sorted.bam",
        "-o", salmon_output_dir
    ]
    logging.info(salmon_command)
    subprocess.run(salmon_command, check=True)
    logging.info(f"Salmon quantification completed. Results saved to: {salmon_output_dir}")

    # Step 2: Move quant.sf to transcript_counts.tsv
    quant_sf_path = os.path.join(salmon_output_dir, "quant.sf")
    transcript_counts_path = os.path.join(salmon_output_dir, "transcript_counts.tsv")

    if os.path.exists(quant_sf_path):
        shutil.move(quant_sf_path, transcript_counts_path)
        logging.info(f"Moved {quant_sf_path} to {transcript_counts_path}")
    else:
        logging.error(f"Error: {quant_sf_path} not found.")


def cpm_gene_count(tmp_dir: str, prefix: str):
    """
    Processes transcript_counts and transcript-to-gene mapping files 
    to compute gene-level counts and CPM (Counts Per Million).
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :param prefix: (str) The prefix to use for output files.
    """
    transcript_counts_file = f"{tmp_dir}/salmon/{prefix}/transcript_counts.tsv"
    tx2gene_file = f"/opt/attrs.tsv"
    # Load the quant.sf file into a DataFrame
    df = pd.read_csv(transcript_counts_file, sep="\t")
    # Load the transcript-to-gene mapping
    tx2gene = pd.read_csv(tx2gene_file, sep="\t")
    # Ensure 'NumReads' is numeric (in case it's in a non-numeric format)
    df['NumReads'] = pd.to_numeric(df['NumReads'], errors='coerce')
    # Merge the quant.sf data with the tx2gene mapping (merge on transcript_id)
    merged = pd.merge(df, tx2gene, left_on="Name", right_on="transcriptId", how="left")
    # Group by geneId and sum the NumReads to get gene counts
    gene_counts = merged.groupby(["geneId", "geneName"])["NumReads"].sum().reset_index()
    # Calculate the total number of reads in the dataset
    total_reads = merged["NumReads"].sum()
    # Calculate CPM (Counts Per Million)
    gene_counts["CPM"] = (gene_counts["NumReads"] / total_reads) * 1e6
    # gene_counts["TPM"] = gene_counts["TPM"].round().astype(int)
    gene_counts.columns = ['geneId', 'geneName', 'NumReads', 'CPM']
    # Save the results as a .tsv file
    os.makedirs("/output", exist_ok=True)
    output_file = f"/output/{prefix}_gene_counts_cpm.tsv"
    gene_counts.to_csv(output_file, sep="\t", index=False)
    logging.info(f"Gene counts and CPM saved to {output_file}")

    # Create the OTTER input file
    gene_results = gene_counts[['geneId', 'geneName', 'NumReads', 'CPM']].copy()

    gene_results = gene_results.rename(columns={
        'geneId': 'gene_id',
        'geneName': 'gene_name',
        'NumReads': 'expected_count',  # This will be overwritten by the new column anyway
        'CPM': 'TPM'
    })
    
    # Save the second results .tsv file
    results_file = f"/output/{prefix}.genes.results"
    gene_results.to_csv(results_file, sep="\t", index=False)
    logging.info(f"OTTER input gene results saved to {results_file}")


def submit_to_otter(token: str, tpm: str, email: str = None, save: bool = False,
                    om:bool = False, v1: bool = False) -> None:
    """
    Submits the data TPM file to the otter web app throug hthe API.
    :param token: (str) API token for the account to be used for submission.
    :param tpm: (str) The filepath to the RSEM TPM file.
    :param email: (list) A list of email addresses that the otter web app will notify of results completion.
    :param save: (bool) A switch to allow the otter web app to save the TPM file sent to it.
    :param om: (bool) A switch to use the otter web app's otter model.
    :param v1: (bool) A switch to use the otter web app's version 1 atlas.
    :return: None
    """
    df = pd.read_csv(tpm, sep='\t')
    data_list = df.to_dict(orient='list')
    data_name = tpm.split('/')[-1]

    headers = {"Authorization": f"Bearer {token}"}
    post_data = {"version": "hierarchical", "data": data_list, "name": data_name, "save": save, "atlas_version": "v2",}
    if om:
        post_data["version"] = "otter"
        post_data["atlas_version"] = "v1"
    if v1:
        post_data["atlas_version"] = "v1"
    if email:
        email_split = email.split(",") if email is not None else None
        post_data["share_with"] = email_split
    for key, value in post_data.items():
        if key != "data":
            logging.debug(f"{key}: {value}")
    r = reqs.post("https://otter.ccm.sickkids.ca/api/inference", json=post_data, headers=headers)
    if r.status_code == 200:
        # Worked
        res = json.loads(r.text)
        this_task_id = res["task_id"]
        success = False
        while not success:
            r_check = reqs.get(f"https://otter.ccm.sickkids.ca/api/inference_check?task_id={this_task_id}",
                               headers=headers)
            if r_check.status_code == 200:
                # Worked
                res_check = json.loads(r_check.text)
                logging.debug(res_check)
                if res_check["position"] is None:
                    # Finished running
                    success = True
                elif res_check["position"] == -1:
                    # Still waiting to be placed in queue
                    pass
                elif res_check["position"] == 0:
                    # Currently running
                    pass
                elif res_check["position"] > 0:
                    # Waiting in queue
                    pass
                else:
                    raise Exception("Problem with position in queue number")
            else:
                logging.info(r_check.status_code)
                logging.debug(r_check.text)
                break
            time.sleep(2)
    else:
        logging.info(f"Failed to post {data_name}")
        logging.info(r.status_code, r.text)
        exit(1)
    return None


def main(config: Any, files: str, tmp_dir: str, prefix: str, cdna: str, drna: bool, 
         token: str, email: str, save: bool, om: bool = False, v1: bool = False) -> None:
    """
    Run the workflow of tasks
    :param files: (str) filepaths of read1 files.
    :param tmp_dir: (str) path of temporary directory for storing intermediate files.
    :param prefix: (str) prefix used on output files.
    :param cdna: (str) cDNA reads library kit.
    :param drna: (bool)  if set to direct RNA sequencing.
    :param token: (str) API token for the otter web app.
    :param email: (str) Email addresses to notify when otter results are ready.
    :param save: (bool) Flag for telling otter web app to save the TPM file data.
    :return: None
    """
    # Get default options
    ref_transcriptome = os.path.join("/opt/", config.get('ref_transcriptome', None))
    fastcat_opts = config.get('fastcat_opts', "-x")
    pychopper_backend = config.get('pychopper_backend', "edlib")
    minimap2_index_opts = config.get('minimap2_index_opts', "-k 14")
    minimap2_cdna_tx_opts = config.get('minimap2_cdna_tx_opts', "-ax map-ont -p 1.0")
    minimap2_drna_tx_opts = config.get('minimap2_drna_tx_opts', "-ax map-ont --for-only -p 1.0")
    aln_filtering_opts = config.get('aln_filtering_opts', "")

    # Validate input files
    inputted_files = [file.strip() for file in files.split(',')]
    try: 
        logging.info("VALIDATING INPUT FILES")
        exact_paths = validate_input_files(inputted_files)
    except ValueError as e:
        logging.error(e)
        exit(1)

    used_fastcat = True
    if len(inputted_files) == 1 and '.fastq' not in inputted_files[0]:
        logging.info("CONCATENATING INPUT FILES USING FASTCAT")
        concatenation(tmp_dir, prefix, inputted_files, fastcat_opts)
    elif len(inputted_files) > 1:
        logging.info("CONCATENATING INPUT FILES USING FASTCAT")
        concatenation(tmp_dir, prefix, exact_paths, fastcat_opts)
    else:
        used_fastcat = False

    reference_indexing(tmp_dir, minimap2_index_opts, ref_transcriptome)
    if drna is False:
        # run pychopper
        logging.info("TRIMMING READS USING PYCHOPPER")
        adaptor_trimming(tmp_dir, used_fastcat, exact_paths[0], cdna, pychopper_backend, prefix)
        # run minimap2
        logging.info("TRANSCRIPTOME MAPPING USING MINIMAP2")
        cdna_transcriptome_mapping(tmp_dir, minimap2_cdna_tx_opts, aln_filtering_opts, prefix)
    else:
        # run minimap2
        logging.info("TRANSCRIPTOME MAPPING USING MINIMAP2")
        drna_transcriptome_mapping(tmp_dir, used_fastcat, exact_paths[0], minimap2_drna_tx_opts, aln_filtering_opts, prefix)

    # run salmon
    logging.info("QUANTIFYING USING SALMON")
    quantitation(tmp_dir, ref_transcriptome, prefix)

    # output .tsv file
    logging.info("OUTPUTTING CPM GENE COUNTS")
    cpm_gene_count(tmp_dir, prefix)

    otter_input_file_path = f"/output/{prefix}.genes.results"
    if os.path.exists(otter_input_file_path) and token is not None:
        logging.info("SUBMITTING TO OTTER.")
        submit_to_otter(token, otter_input_file_path, email, save, om, v1)


if __name__ == "__main__":
    # Load defaults from config file
    config = load_config('/opt/params.config')

    parser = argparse.ArgumentParser(prog='long-read',
                                     description="Process RNA-seq fastqs to gene expression counts, TPM.")
    parser.add_argument("-i", "--input", required=True, metavar='', type=str,
                        help="input file(s) seperated by comma or single input directory.")
    parser.add_argument("-p", "--prefix", required=True, metavar='', type=str,
                        help="output prefix.")
    parser.add_argument("-c", "--cdna", default=config.get('cdna', "PCB114"), metavar='', type=str,
                        help="cDNA reads library kit used, SQK-PCB114 by default.")
    parser.add_argument("-d", "--drna", action='store_true', dest='drna', default=config.getboolean('drna', False),
                        help="direct RNA sequencing skipping read trimming, set to FALSE by default.")
    parser.add_argument("-t", "--token", metavar='', type=str, help="Otter web app API token.")
    parser.add_argument("-e", "--email", type=str, metavar='',
                        help="Comma separate list of email addresses, otter web app will email when analysis complete.")
    parser.add_argument("-s", "--save", action='store_true', dest='save', default=False,
                        help="NOTE: Setting this will allow the otter web app to keep the TPM file sent.")
    args = parser.parse_args()

    token = args.token
    email = args.email

    prefix = args.prefix
    if prefix is None:
        prefix = "sample"

    cdna = args.cdna
    drna = args.drna
    save = args.save

    # Check /tmp is available
    with tempfile.TemporaryDirectory(prefix="moose_") as tmp_dir_name:
        logging.info(f"Created temporary directory: {tmp_dir_name}")
        main(config, tmp_dir_name, prefix, cdna, drna, token, email, save)
    exit(0)
