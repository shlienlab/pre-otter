#!/usr/bin/env python3

import argparse
import logging
import configparser
import sys
import os
import pandas as pd
import requests as reqs
import tempfile

from scripts import short_read, long_read

# Set up logging configuration to log to both terminal and file
logging_dir = "/output"
if not os.path.exists(logging_dir):
    logging_dir = os.getcwd()
    print(f"/output directory not found for logging. Using current directory instead: {logging_dir}")

try:
    logging.basicConfig(
        level=logging.INFO,  # Minimum level to capture
        format='%(asctime)s - %(levelname)s - %(message)s',  # Log format
        datefmt='%Y-%m-%d %H:%M:%S %p %Z',  # Date format
        handlers=[
            logging.FileHandler(os.path.join(logging_dir, "pre-otter.log")),   # Log to a file
            logging.StreamHandler()           # Log to the console (terminal)
        ]
    )
    logger = logging.getLogger()
except FileNotFoundError as e:
    print(f"Error configuring logging: {e}")
    sys.exit(1)


def set_log_level(debug: bool):
    """Function to set the logging level dynamically."""
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

# Common argument for both modes
def common_args(parser):
    parser.add_argument("-p", "--prefix", required=True, metavar='', type=str,
                        help="Prefix to be used for the output files.")
    parser.add_argument("--om", action="store_true", required=False,
                        help="Use Otter model instead of hierarchical model.",default=False)
    parser.add_argument("--v1", required=False, action="store_true",
                        help="Use version 1 of model instead of version 2.",default=False)
    parser.add_argument("-t", "--token", metavar='', type=str, help="Otter web app API token.")
    parser.add_argument("-e", "--email", type=str, metavar='',
                        help="Comma separate list of email addresses, otter web app will email when analysis complete.")
    parser.add_argument("-s", "--save", action='store_true', dest='save', default=False,
                        help="NOTE: Setting this will allow the otter web app to keep the TPM file sent.")


def load_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    return config['defaults']  # Return the default values from the config file


def parse_filenames(value):
    # This regular expression will match a list of filenames, separated by commas and possibly spaces
    files = value.split(',')
    return [item for item in files if item != '']


def run_long_read(config, args):
    with tempfile.TemporaryDirectory(prefix="otter_") as tmp_dir_name:
        logging.info(f"Created temporary directory: {tmp_dir_name}")
        long_read.main(config=config, files=args.input, tmp_dir=tmp_dir_name, prefix=args.prefix,
                       cdna=args.cdna, drna=args.drna, token=args.token, email=args.email, 
                       save=args.save, om=args.om, v1=args.v1)
    logging.info("Returned from long-read script.")


def run_short_read(args):
    with tempfile.TemporaryDirectory(prefix="otter_") as tmp_dir_name:
        logging.info(f"Created temporary directory: {tmp_dir_name}")
        short_read.main(this_read1s=args.read1, this_read2s=args.read2, tmp_dir=tmp_dir_name, prefix=args.prefix,
                        token=args.token, email=args.email, save=args.save, om=args.om, v1=args.v1)
    logging.info("Returned from short-read script.")


def main():
    # Load defaults from config file
    try:
        config = load_config('/opt/params.config')
    except KeyError:
        config = load_config(os.path.join(repo_dir, "scripts",  "params.config"))
    parser = argparse.ArgumentParser(prog='pre-otter',
                                     description="Process RNA-seq fastqs to gene expression counts, TPM.")

    parser.add_argument("-b", "--debug", action='store_true', dest='debug', default=False, help="Debug.")
    # General mode argument to select long-read or short-read
    # Define the subparsers for different modes
    subparsers = parser.add_subparsers(dest='mode', help='Mode of operation')
    short_read_parser = subparsers.add_parser('short-read', help='Run in short-read mode.')
    long_read_parser = subparsers.add_parser('long-read', help='Run in long-read mode.')

    # Arguments for short-read mode
    short_read_parser.add_argument("-r1", "--read1", required=True, metavar='', type=str,
                        help="fastq read1 files, comma separated if more than one.")
    short_read_parser.add_argument("-r2", "--read2", required=True, metavar='', type=str,
                        help="fastq read2 files, comma separated if more than one.")
    common_args(short_read_parser)

    # Arguments for long-read mode
    long_read_parser.add_argument("-i", "--input", required=True, metavar='', type=str,
                        help="input file(s) seperated by comma or single input directory.")
    long_read_parser.add_argument("-c", "--cdna", default=config.get('cdna', "PCB114"), metavar='', type=str,
                        help="cDNA reads library kit used, SQK-PCB114 by default.")
    long_read_parser.add_argument("-d", "--drna", action='store_true', dest='drna', default=config.getboolean('drna', False),
                        help="direct RNA sequencing skipping read trimming, set to FALSE by default.")
    common_args(long_read_parser)

    # Parse all arguments
    args = parser.parse_args()

    if args.debug:
        set_log_level(debug=True)

    # Handle the case where mode is not specified
    if args.mode is None:
        print("Please specify a mode: short-read or long-read.")
        return
    
    # Warn about lack of user credentials
    if args.token is None:
        logging.warning("OTTER token was not specified.")

    if args.email is None:
        logging.warning("User email was not specified.")

    # Handle mode selection and script execution
    if args.mode == "long-read":
        if args.input is None:
            print("For long-read mode, -i is required.")
            sys.exit(1)
        run_long_read(config, args)

    elif args.mode == "short-read":
        if args.read1 is None or args.read2 is None:
            print("For short-read mode, -r1 and -r2 are required.")
            sys.exit(1)
        run_short_read(args)


if __name__ == "__main__":
    repo_dir = os.path.dirname(os.path.realpath(__file__))
    main()
    exit(0)
