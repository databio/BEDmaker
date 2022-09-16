#!/usr/bin/env python3

from argparse import ArgumentParser
import pypiper
import os
# import re
import sys
import tempfile
import pandas as pd
import gzip
import shutil
from refgenconf import (
    RefGenConf as RGC,
    select_genome_config,
    RefgenconfError,
    CFG_ENV_VARS,
    CFG_FOLDER_KEY,
)
from yacman.exceptions import UndefinedAliasError

# CONSTANTS
# Creating list of standard chromosome names:
STANDARD_CHROM_LIST = ["chr" + str(chr_nb) for chr_nb in list(range(1, 23))]
STANDARD_CHROM_LIST[len(STANDARD_CHROM_LIST):] = ["chrX", "chrY", "chrM"]

parser = ArgumentParser(description="A pipeline to convert bigwig or bedgraph files into bed format")

parser.add_argument("-f", "--input-file", help="path to the input file", type=str)
parser.add_argument(
    "-n",
    "--narrowpeak",
    help="whether the regions are narrow (transcription factor implies narrow, histone mark implies broad peaks)",
    type=bool,
)
parser.add_argument(
    "-t",
    "--input-type",
    help="a bigwig or a bedgraph file that will be converted into BED format",
    type=str,
)
parser.add_argument("-g", "--genome", help="reference genome", type=str)
parser.add_argument("-r", "--rfg-config", help="file path to the genome config file", type=str)
parser.add_argument("-o", "--output-bed", help="path to the output BED files", type=str)
parser.add_argument("--output-bigbed", help="path to the output bigBed files", type=str)
parser.add_argument(
    "-s",
    "--sample-name",
    help="name of the sample used to systematically build the output name",
    type=str,
)
parser.add_argument(
    "--chrom-sizes",
    help="a full path to the chrom.sizes required for the bedtobigbed conversion",
    type=str,
    required=False
)
parser.add_argument(
    "--standard-chrom",
    help="Standardize chromosome names. Default: False",
    action="store_true"
)
# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(
    parser,
    groups=["pypiper", "looper"],
    required=["--input-file", "--input-type"]
)

args = parser.parse_args()

# COMMANDS TEMPLATES

# bedGraph to bed
bedGraph_template = "macs2 {width} -i {input} -o {output}"
# bigBed to bed
bigBed_template = "bigBedToBed {input} {output}"
# bigWig to bed
bigWig_template = (
    "bigWigToBedGraph {input} /dev/stdout | macs2 {width} -i /dev/stdin -o {output}"
)
# preliminary for wig to bed
# wig_template =  "wigToBigWig {input} {chrom_sizes} /dev/stdout -clip | bigWigToBedGraph /dev/stdin  /dev/stdout | macs2 {width} -i /dev/stdin -o {output}"
wig_template = "wigToBigWig {input} {chrom_sizes} {intermediate_bw} -clip"
# bed default link
# bed_template = "ln -s {input} {output}"
bed_template = "cp {input} {output}"
# gzip output files
gzip_template = "gzip {unzipped_converted_file} "

# SET OUTPUT FOLDERS

file_name = os.path.basename(args.input_file)
file_id = os.path.splitext(file_name)[0]
input_extension = os.path.splitext(file_name)[1]  # is it gzipped or not?
output_bed_extension = os.path.splitext(args.output_bed)[1]
# sample_folder = os.path.join(out_parent, args.sample_name)  # specific output folder for each sample log and stats
if args.input_type != "bed":
    if input_extension == ".gz":
        output_bed = (
                os.path.splitext(os.path.splitext(args.output_bed)[0])[0] + ".bed.gz"
        )
    else:
        output_bed = os.path.splitext(args.output_bed)[0] + ".bed.gz"
else:
    if input_extension != ".gz" and output_bed_extension != ".gz":
        output_bed = args.output_bed + ".gz"
    else:
        output_bed = args.output_bed

bed_parent = os.path.dirname(args.output_bed)
if not os.path.exists(bed_parent):
    print(f"Output directory does not exist. Creating: {bed_parent}")
    os.makedirs(bed_parent)

if not os.path.exists(args.output_bigbed):
    print(f"BigBed directory does not exist. Creating: {args.output_bigbed}")
    os.makedirs(args.output_bigbed)

logs_name = "bedmaker_logs"
logs_dir = os.path.join(bed_parent, logs_name, args.sample_name)
if not os.path.exists(logs_dir):
    print("bedmaker logs directory doesn't exist. Creating one...")
    os.makedirs(logs_dir)

pm = pypiper.PipelineManager(
    name="bedmaker", outfolder=logs_dir, args=args
)  # ArgParser and add_pypiper_args


def get_chrom_sizes():
    """
    Get chrom.sizes file path by input arg, or Refegenie.

    :return str: chrom.sizes file path
    """

    if args.chrom_sizes:
        return args.chrom_sizes

    print("Determining path to chrom.sizes asset via Refgenie.")
    # get path to the genome config; from arg or env var if arg not provided
    refgenie_cfg_path = select_genome_config(
        filename=args.rfg_config, check_exist=False
    )
    if not refgenie_cfg_path:
        raise OSError(
            "Could not determine path to a refgenie genome configuration file. "
            "Use --rfg-config argument or set '{}' environment variable to provide it".format(
                CFG_ENV_VARS
            )
        )
    if isinstance(refgenie_cfg_path, str) and not os.path.exists(refgenie_cfg_path):
        # file path not found, initialize a new config file
        print(f"File '{refgenie_cfg_path}' does not exist. Initializing refgenie genome configuration file.")
        rgc = RGC(entries={CFG_FOLDER_KEY: os.path.dirname(refgenie_cfg_path)})
        rgc.initialize_config_file(filepath=refgenie_cfg_path)
    else:
        print(f"Reading refgenie genome configuration file from file: {refgenie_cfg_path}")
        rgc = RGC(filepath=refgenie_cfg_path)
    try:
        # get path to the chrom.sizes asset
        chrom_sizes = rgc.seek(
            genome_name=args.genome,
            asset_name="fasta",
            tag_name="default",
            seek_key="chrom_sizes",
        )
        print(chrom_sizes)
    except (UndefinedAliasError, RefgenconfError):
        # if chrom.sizes not found, pull it first
        print("Could not determine path to chrom.sizes asset, pulling")
        rgc.pull(genome=args.genome, asset="fasta", tag="default")
        chrom_sizes = rgc.seek(
            genome_name=args.genome,
            asset_name="fasta",
            tag_name="default",
            seek_key="chrom_sizes",
        )

    print("Determined path to chrom.sizes asset: {}".format(chrom_sizes))

    return chrom_sizes


def get_bed_type(bed):
    """
    get bed type

    :return bed type
    """
    #    column format for bed12
    #    string chrom;       "Reference sequence chromosome or scaffold"
    #    uint   chromStart;  "Start position in chromosome"
    #    uint   chromEnd;    "End position in chromosome"
    #    string name;        "Name of item."
    #    uint score;          "Score (0-1000)"
    #    char[1] strand;     "+ or - for strand"
    #    uint thickStart;   "Start of where display should be thick (start codon)"
    #    uint thickEnd;     "End of where display should be thick (stop codon)"
    #    uint reserved;     "Used as itemRgb as of 2004-11-22"
    #    int blockCount;    "Number of blocks"
    #    int[blockCount] blockSizes; "Comma separated list of block sizes"
    #    int[blockCount] chromStarts; "Start positions relative to chromStart"
    df = pd.read_csv(bed, sep="\t", header=None)
    df = df.dropna(axis=1)

    # standardizing chromosome names
    if args.standard_chrom:
        print("Standardizing chromosomes...")
        df = df[df.loc[:, 0].isin(STANDARD_CHROM_LIST)]
        df.to_csv(bed, compression='gzip', sep="\t", header=False, index=False)

    num_cols = len(df.columns)
    print(df)
    bedtype = 0
    for col in df:
        if col <= 2:
            if col == 0:
                if df[col].dtype == "O":
                    bedtype += 1
                else:
                    return None
            else:
                if df[col].dtype == "int" and (df[col] >= 0).all():
                    bedtype += 1
                else:
                    return None
        else:
            if col == 3:
                if df[col].dtype == "O":
                    bedtype += 1
                else:
                    n = num_cols - bedtype
                    return f"bed{bedtype}+{n}"
            elif col == 4:
                if df[col].dtype == "int" and df[col].between(0, 1000).all():
                    bedtype += 1
                else:
                    n = num_cols - bedtype
                    return f"bed{bedtype}+{n}"
            elif col == 5:
                if df[col].isin(["+", "-", "."]).all():
                    bedtype += 1
                else:
                    n = num_cols - bedtype
                    return f"bed{bedtype}+{n}"
            elif 6 <= col <= 8:
                if df[col].dtype == "int" and (df[col] >= 0).all():
                    bedtype += 1
                else:
                    n = num_cols - bedtype
                    return f"bed{bedtype}+{n}"
            elif col == 9:
                if df[col].dtype == "int":
                    bedtype += 1
                else:
                    n = num_cols - bedtype
                    return f"bed{bedtype}+{n}"
            elif col == 10 or col == 11:
                if df[col].str.match("^(\d+(,\d+)*)?$").all():
                    bedtype += 1
                else:
                    n = num_cols - bedtype
                    return f"bed{bedtype}+{n}"
            else:
                n = num_cols - bedtype
                return f"bed{bedtype}+{n}"


def main():
    # pm = pypiper.PipelineManager(name="bedmaker", outfolder=logs_dir, args=args) # ArgParser and add_pypiper_args

    # Define target folder for converted files and implement the conversions; True=TF_Chipseq False=Histone_Chipseq

    print(f"Got input type: {args.input_type}")
    if not args.input_type == "bed":
        print(f"Converting {args.input_file} to BED format")

    # Define whether Chip-seq data has broad or narrow peaks
    width = "bdgbroadcall" if not args.narrowpeak else "bdgpeakcall"

    # # Call pyBigWig to ensure bigWig and bigBed files have the correct format
    # if args.input_type in ["bigWig", "bigBed"]:
    #     obj = pyBigWig.open(args.input_file)
    #     validation_method = getattr(obj, "isBigBed" if args.input_type == "bigBed" else "isBigWig")
    #     if not validation_method():
    #         raise Exception("{} file did not pass the {} format validation".
    #                         format(args.input_file, args.input_type))

    input_file = args.input_file
    # Use the gzip and shutil modules to produce temporary unzipped files
    if input_extension == ".gz":
        input_file = os.path.join(
            os.path.dirname(args.output_bed), os.path.splitext(file_name)[0]
        )
        with gzip.open(args.input_file, "rb") as f_in:
            with open(input_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        pm.clean_add(input_file)

    temp_bed_path = os.path.splitext(output_bed)[0]
    if args.input_type == "bedGraph":
        cmd = bedGraph_template.format(
            input=input_file, output=temp_bed_path, width=width
        )
    elif args.input_type == "bigWig":
        cmd = bigWig_template.format(
            input=input_file, output=temp_bed_path, width=width
        )
    elif args.input_type == "wig":

        chrom_sizes = get_chrom_sizes()

        # define a target for temporary bw files
        temp_target = os.path.join(bed_parent, file_id + ".bw")
        pm.clean_add(temp_target)
        cmd1 = wig_template.format(
            input=input_file,
            intermediate_bw=temp_target,
            chrom_sizes=chrom_sizes,
            width=width,
        )
        cmd2 = bigWig_template.format(
            input=temp_target, output=temp_bed_path, width=width
        )
        cmd = [cmd1, cmd2]
    elif args.input_type == "bigBed":
        cmd = bigBed_template.format(input=input_file, output=temp_bed_path)
    elif args.input_type == "bed":
        if input_extension == ".gz":
            cmd = bed_template.format(input=args.input_file, output=output_bed)
        else:
            # cmd = [gzip_template.format(unzipped_converted_file=input_file),
            #        bed_template.format(input=input_file + ".gz", output=output_bed)]
            cmd = [bed_template.format(input=input_file, output=os.path.splitext(output_bed)[0]),
                   gzip_template.format(unzipped_converted_file=os.path.splitext(output_bed)[0])]

    else:
        raise NotImplementedError(f"'{args.input_type}' format is not supported")

    if args.input_type != "bed" and input_extension != ".gz":
        gzip_cmd = gzip_template.format(unzipped_converted_file=temp_bed_path)
        if not isinstance(cmd, list):
            cmd = [cmd]
        cmd.append(gzip_cmd)
    pm.run(cmd, target=output_bed)

    print(f"Generating bigBed files for {args.input_file}")

    bedfile_name = os.path.split(output_bed)[1]
    fileid = os.path.splitext(os.path.splitext(bedfile_name)[0])[0]
    # Produce bigBed (big_narrow_peak) file from peak file
    big_narrow_peak = os.path.join(args.output_bigbed, fileid + ".bigBed")

    chrom_sizes = get_chrom_sizes()

    temp = os.path.join(args.output_bigbed, next(tempfile._get_candidate_names()))

    if not os.path.exists(big_narrow_peak):
        bedtype = get_bed_type(output_bed)
        pm.clean_add(temp)
        
        if bedtype is not None:
            cmd = "zcat " + output_bed + "  | sort -k1,1 -k2,2n > " + temp
            pm.run(cmd, temp)
            cmd = f"bedToBigBed -type={bedtype} {temp} {chrom_sizes} {big_narrow_peak}"
            try:
                pm.run(cmd, big_narrow_peak, nofail=True)
            except Exception as err:
                print(
                    f"Fail to generating bigBed files for {args.input_file}: "
                    f"unable to validate genome assembly with Refgenie. Error: {err}"
                )
        else:
            cmd = "zcat " + output_bed + " | awk '{ print $1, $2, $3 }'| sort -k1,1 -k2,2n > " + temp
            pm.run(cmd, temp)
            cmd = f"bedToBigBed -type=bed3 {temp} {chrom_sizes} {big_narrow_peak}"
            try:
                pm.run(cmd, big_narrow_peak, nofail=True)
            except Exception as err:
                print(
                    f"Fail to generating bigBed files for {args.input_file}: "
                    f"unable to validate genome assembly with Refgenie. Error: {err}"
                )
            # print(f"Fail to generating bigBed files for {args.input_file}: invalid bed format")

    pm.stop_pipeline()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
