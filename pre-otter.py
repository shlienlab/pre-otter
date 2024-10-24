#!/usr/bin/env python3

import argparse
import asyncio
import gzip
import multiprocessing
import os
from typing import Tuple

import pandas as pd
import requests as reqs
import subprocess
import shutil
import tempfile

from pathlib import Path

REF_GENOME_VERSION = "hg38"

def check_files_exist(fastqs: str) -> list[str] | None:
    """
    Test that fastq files exist.
    :param fastqs: (str) Filepath to fastq file.
    :return: list of file path strings.
    """
    files = fastqs.split(',')
    for file in files:
        if not os.path.exists(file):
            return None
    return files


async def uncompress_file(file: str, out_dir: str) -> str:
    """
    Used to uncompress gz files required by the pipeline held in the container.
    :param file: (str) File path to the file to be uncompressed.
    :param out_dir: (str) directory for the file to be uncompressed to.
    :return: (str) The uncompressed file.
    """
    this_file_base_no_ext = Path(file).stem
    output_file = os.path.join(out_dir, this_file_base_no_ext)
    if os.path.isfile(output_file):
        # file already exists, return now
        return output_file
    with gzip.open(file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return output_file


async def concat_fastqs(files: list, tmp_dir: str, output_file: str) -> str:
    """
    Concatenate fastq files into a single file.
    :param files: (list) file paths of fqstq files.
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :param output_file: (str) Name of the output file.
    :return: (str) Path to the output file.
    """
    output_file = os.path.join(tmp_dir, output_file)
    with open(output_file, 'wb') as outfile:
        for file in files:
            with open(file, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
    return output_file


def run_fastp(this_fastq1: str, this_fastq2: str, tmp_dir: str) -> tuple[str, str] | None:
    """
    Run the fastp tool to trim the reads
    :param this_fastq1: (str) filepath of read1 file.
    :param this_fastq2: (str) filepath of read2 file.
    :param tmp_dir: (str) Temporary directory for storing intermediate files.
    :return: (tuple) Read1 and Read2 file that have been trimmed.
    """
    out_fq1 = f"{tmp_dir}/fastp_out_R1.fastq.gz"
    out_fq2 = f"{tmp_dir}/fastp_out_R2.fastq.gz"

    command = [f"fastp -l 35 --in1 {this_fastq1} --out1 {out_fq1} --in2 {this_fastq2} --out2 {out_fq2}"]
    print(command)
    completion = subprocess.run(command, shell=True, check=True)
    if completion.returncode != 0:
        return None
    return out_fq1, out_fq2


async def generate_genome(genome_dir: str, fasta: str, gtf: str) -> str | None:
    """
    Run the start command to generate the STAR index used for alignment
    :param genome_dir: (str) The directory to store the index
    :param fasta: (str) File path of the reference fasta file.
    :param gtf: (str) File path of the reference GTF file.
    :return: (str) The directory path of the index
    """
    # First check to see if the files already exist.
    file_missing = False
    files_required = ["chrLength.txt", "chrName.txt", "exonGeTrInfo.tab", "geneInfo.tab", "genomeParameters.txt", "SA",
                      "sjdbInfo.txt", "sjdbList.out.tab", "chrNameLength.txt", "chrStart.txt", "exonInfo.tab", "Genome",
                      "SAindex", "sjdbList.fromGTF.out.tab", "transcriptInfo.tab"]
    for file in files_required:
        if not os.path.isfile(os.path.join(genome_dir, file)):
            file_missing = True
    if not file_missing:
        return genome_dir

    # They don't so carry on
    cpu_count = multiprocessing.cpu_count()

    command = [f"STAR --runThreadN {cpu_count} --runMode genomeGenerate --genomeDir {genome_dir} "
               f"--genomeFastaFiles {fasta} --sjdbGTFfile {gtf}"]
    print(command)
    completion = subprocess.run(command, shell=True, check=True)
    if completion.returncode != 0:
        return None
    print(next(os.walk(genome_dir))[2])
    return genome_dir


# def decompress_gz_files(root_dir):
#     """
#
#     :param root_dir:
#     :return:
#     """
#     for root, dirs, files in os.walk(root_dir):
#         for file in files:
#             if file.endswith('.gz'):
#                 file_path = os.path.join(root, file)
#                 with gzip.open(file_path, 'rb') as f_in:
#                     with open(file_path[:-3], 'wb') as f_out:  # Decompressed file will have the original name without the .gz extension
#                         f_out.write(f_in.read())
#                 os.remove(file_path)  # Remove the original compressed file


async def rsem_prepare_reference(ref_dir: str, fasta: str, gtf: str):
    """
    Prepare the files that RSEM requires to run.
    :param ref_dir: (str) The directory to store the files.
    :param fasta: (str) The path to the reference fasta file.
    :param gtf: (str) The path to the reference GTF file.
    :return: (str) The path ot the directory the RSEM files are stored.
    """
    ## First check to see if the files already exist.
    file_missing = False
    files_required = [f"{REF_GENOME_VERSION}.chrlist", f"{REF_GENOME_VERSION}.grp", f"{REF_GENOME_VERSION}.idx.fa",
                      f"{REF_GENOME_VERSION}.n2g.idx.fa", f"{REF_GENOME_VERSION}.seq", f"{REF_GENOME_VERSION}.ti",
                      f"{REF_GENOME_VERSION}.transcripts.fa"]
    for file in files_required:
        if not os.path.isfile(os.path.join(ref_dir, file)):
            file_missing = True
    if not file_missing:
        return ref_dir

    # They don't so carry on

    cpu_count = multiprocessing.cpu_count()
    command = [f"rsem-prepare-reference --gtf {gtf} --num-threads {cpu_count} {fasta} {ref_dir}/{REF_GENOME_VERSION}"]
    print(command)
    completion = subprocess.run(command, shell=True, check=True)
    if completion.returncode != 0:
        return None
    print(next(os.walk(ref_dir))[2])
    return ref_dir


def run_star(genome_dir: str, read1: str, read2: str, tmp_dir: str):
    """
    Run the star aligner on the fastq files.
    :param genome_dir: (str) The directory to store the index
    :param read1: (str) filepath of read1 file.
    :param read2: (str) filepath of read2 file.
    :param tmp_dir: (str) The temporary directory for intermediate files.
    :return: (str) The path ot the bam file output of STAR.
    """
    cpu_count = multiprocessing.cpu_count()
    command = [f"STAR --runMode alignReads --runThreadN {cpu_count} --genomeDir {genome_dir} "
               f"--readFilesIn {read1} {read2} --outFileNamePrefix {tmp_dir}/moose_rna --outSAMunmapped Within "
               f"--quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD --outFilterType BySJout "
               f"--outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 "
               f"--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 "
               f"--alignSJDBoverhangMin 1 --sjdbScore 1 --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
               f"--limitOutSJcollapsed 5500000 --limitBAMsortRAM 48000000000"]
    print(command)
    completion = subprocess.run(command, shell=True, check=True)
    if completion.returncode != 0:
        return None
    return f"{tmp_dir}/moose_rnaAligned.toTranscriptome.out.bam"


def run_rsem(tmp_dir: str, bam_file: str, output_dir: str, prefix: str):
    """
    Run the RSEM expression calculations.
    :param tmp_dir: (str) The temporary directory to store intermediate files.
    :param bam_file: (str) The filepath to the bam file to have reads counted.
    :param output_dir: (str) The output directory to store results.
    :param prefix: (str) The prefix to use for output files created by RSEM.
    :return: (str) The path to the *.genes.results file from RSEM.
    """
    cpu_count = multiprocessing.cpu_count()
    command = [f"rsem-calculate-expression --quiet --no-qualities --num-threads {cpu_count} --paired-end "
               f"--forward-prob 0.5 --seed-length 25 --fragment-length-mean -1 --no-bam-output --bam {bam_file} "
               f"{tmp_dir}/{REF_GENOME_VERSION} /output/{prefix}"]
    print(command)
    completion = subprocess.run(command, shell=True, check=True)
    if completion.returncode != 0:
        return None
    return f"/output/{prefix}.genes.results"


def submit_to_otter(token: str, tpm: str, email: list = None, save: bool = False) -> None:
    """
    Submits the data TPM file to the otter web app throug hthe API.
    :param token: (str) API token for the account to be used for submission.
    :param tpm: (str) The filepath to the RSEM TPM file.
    :param email: (list) A list of email addresses that the otter web app will notify of results completion.
    :param save: (bool) A switch to allow the otter web app to save the TPM file sent to it.
    :return: None
    """
    df = pd.read_csv(tpm, sep='\t')
    data = df.to_dict(orient='list')
    name = tpm.split('/')[-1]

    headers = {'Authorization': f'Bearer {token}'}
    if email is None:
        post_data = {'version': 'otter', 'data': data, 'name': name, "save": save}
    else:
        post_data = {'version': 'otter', 'data': data, 'name': name, 'share_with': email, "save": save}
    r = reqs.post('https://otter.ccm.sickkids.ca/api/inference', json=post_data, headers=headers)
    if r.status_code == 200:
        print(f"{name} have been submitted to otter.")
    else:
        print(r.status_code, r.text)
    return None


async def main(this_read1s: list, this_read2s: list, tmp_dir: str, prefix: str, token: str, save: bool, email: list) -> None:
    """
    Run the workflow of tasks
    :param this_read1s: (list) filepaths of read1 files.
    :param this_read2s: (list) filepaths of read2 files.
    :param tmp_dir: (str) path of temporary directory for storing intermediate files.
    :param prefix: (str) prefix used on output files.
    :param token: (str) API token for the otter web app.
    :param save: (bool) Flag for telling otter web app to save the TPM file data.
    :param email: (list) Email addresses to notify when otter results are ready.
    :return: None
    """
    gtf = lambda: uncompress_file(file="/opt/gencode.v23.annotation.gtf.gz", out_dir="/reference")
    fasta = lambda: uncompress_file(file="/opt/hg38.fa.gz", out_dir="/reference")

    uncompress_results = await asyncio.gather(
        fasta(),
        gtf(),
    )
    my_fasta = uncompress_results[0]
    my_gtf = uncompress_results[1]

    read1 = lambda: concat_fastqs(files=this_read1s, tmp_dir=tmp_dir, output_file="read1.fastq.gz")
    read2 = lambda: concat_fastqs(files=this_read2s, tmp_dir=tmp_dir, output_file="read2.fastq.gz")

    genome_dir = lambda: generate_genome(genome_dir="/reference", fasta=my_fasta, gtf=my_gtf)
    rsem_dir = lambda: rsem_prepare_reference(ref_dir="/reference", fasta=my_fasta, gtf=my_gtf)

    results = await asyncio.gather(
        read1(),
        read2(),
        genome_dir(),
        rsem_dir(),
    )
    my_read1 = results[0]
    my_read2 = results[1]
    my_references = results[2]
    my_genome_dir = results[3]

    if my_references is None:
        print(f"References and indexes not present.")
        exit(1)

    print(results)

    my_cut_fq1, my_cut_fq2 = run_fastp(this_fastq1=my_read1, this_fastq2=my_read2, tmp_dir=tmp_dir)

    # Now run Star
    bam_file = run_star(genome_dir=my_genome_dir, read1=my_cut_fq1, read2=my_cut_fq2, tmp_dir=tmp_dir)

    # Now run rsem
    if bam_file is not None:
        rsem_out = run_rsem(tmp_dir=my_references, bam_file=bam_file, output_dir="/output", prefix=prefix)
        if rsem_out is not None and token is not None:
            submit_to_otter(token=my_token, tpm=f"/output/{prefix}.genes.results", email=email, save=save)

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='pre-otter',
                                     description="Process RNA-seq fastqs to gene expression counts, TPM.")
    parser.add_argument("-r1", "--read1", required=True, metavar='', type=str,
                        help="fastq read1 files, comma separated if more than one.")
    parser.add_argument("-r2", "--read2", required=True, metavar='', type=str,
                        help="fastq read2 files, comma separated if more than one.")
    parser.add_argument("-p", "--prefix", required=True, metavar='', type=str,
                        help="Prefix to be used for the output files.")
    parser.add_argument("-t", "--token", metavar='', type=str, help="Otter web app API token.")
    parser.add_argument("-e", "--email", type=str, metavar='',
                        help="Comma separate list of email addresses, otter web app will email when analysis complete.")
    parser.add_argument("-s", "--save", action='store_true', dest='save', default=False,
                        help="NOTE: Setting this will allow the otter web app to keep the TPM file sent.")
    args = parser.parse_args()

    read1s, read2s, my_email = None, None, None
    if args.read1 is not None:
        read1s = check_files_exist(fastqs=args.read1)
    if args.read2 is not None:
        read2s = check_files_exist(fastqs=args.read2)
    my_prefix = args.prefix
    my_token = args.token
    if args.email is not None:
        my_email = args.email.split(",")
    my_save = args.save

    # Check files exist
    if read1s is None or read2s is None:
        print("Files not found.")
        exit(1)

    # Check /tmp is available
    with tempfile.TemporaryDirectory(prefix="moose_") as tmp_dir_name:
        print('created temporary directory', tmp_dir_name)
        asyncio.run(main(this_read1s=read1s, this_read2s=read2s,
                         tmp_dir=tmp_dir_name, prefix=my_prefix,
                         token=my_token, email=my_email, save=my_save))
    exit(0)
