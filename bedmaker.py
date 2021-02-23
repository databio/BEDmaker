#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os
import sys
import tempfile
# import pyBigWig
import gzip
import shutil
from refgenconf import RefGenConf as RGC, select_genome_config, RefgenconfError, CFG_ENV_VARS, CFG_FOLDER_KEY

parser = ArgumentParser(description="A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("-f", "--input-file", help="path to the input file", type=str)
parser.add_argument("-n", "--narrowpeak", help="whether the regions are narrow (transcription factor implies narrow, histone mark implies broad peaks)", type=bool)
parser.add_argument("-t", "--input-type", help="a bigwig or a bedgraph file that will be converted into BED format", type=str)
parser.add_argument("-g", "--genome", help="reference genome", type=str)
parser.add_argument("-r", "--rfg-config", help="file path to the genome config file", type=str)
parser.add_argument("-o", "--output-bed", help="path to the output BED files", type=str)
parser.add_argument("--output-bigbed", help="path to the output bigBed files", type=str)
parser.add_argument("-s", "--sample-name", help="name of the sample used to systematically build the output name", type=str)
parser.add_argument("--chrom-sizes", help="a full path to the chrom.sizes required for the bedtobigbed conversion", type=str)

# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                  required=["--input-file", "--input-type"])

args = parser.parse_args()

# COMMANDS TEMPLATES

# bedGraph to bed
bedGraph_template = "macs2 {width} -i {input} -o {output}"
# bigBed to bed
bigBed_template = "bigBedToBed {input} {output}"
# bigWig to bed
bigWig_template = "bigWigToBedGraph {input} /dev/stdout | macs2 {width} -i /dev/stdin -o {output}"
# preliminary for wig to bed
# wig_template =  "wigToBigWig {input} {chrom_sizes} /dev/stdout -clip | bigWigToBedGraph /dev/stdin  /dev/stdout | macs2 {width} -i /dev/stdin -o {output}"
wig_template = "wigToBigWig {input} {chrom_sizes} {intermediate_bw} -clip" 
# bed default link
bed_template = "ln -s {input} {output}"
# gzip output files
gzip_template = "gzip {unzipped_converted_file} "

# SET OUTPUT FOLDERS

file_name = os.path.basename(args.input_file)
file_id = os.path.splitext(file_name)[0]
input_extension = os.path.splitext(file_name)[1]  # is it gzipped or not?
#sample_folder = os.path.join(out_parent, args.sample_name)  # specific output folder for each sample log and stats

bed_parent = os.path.dirname(args.output_bed)
if not os.path.exists(bed_parent):
    print("Output directory does not exist. Creating: {}".format(bed_parent))
    os.makedirs(bed_parent)

if not os.path.exists(args.output_bigbed):
    print("BigBed directory does not exist. Creating: {}".format(args.output_bigbed))
    os.makedirs(args.output_bigbed)

logs_name = "bedmaker_logs"
logs_dir = os.path.join(bed_parent, logs_name, args.sample_name)
if not os.path.exists(logs_dir):
    print("bedmaker logs directory doesn't exist. Creating one...")
    os.makedirs(logs_dir)

def get_chrom_sizes():
    """
    Get chrom.sizes file path by input arg, or Refegenie.

    :return str: chrom.sizes file path
    """
    
    if args.chrom_sizes:
        chrom_sizes = args.chrom_sizes
    elif args.rfg_config:
        print("Chrom.sizes file is not specified.  Determining path to chrom.sizes asset via Refgenie.")
        # get path to the genome config; from arg or env var if arg not provided
        refgenie_cfg_path = select_genome_config(filename=args.rfg_config, check_exist=False)
        if not refgenie_cfg_path:
            raise OSError("Could not determine path to a refgenie genome configuration file. "
                          "Use --rfg-config argument or set '{}' environment variable to provide it".
                          format(CFG_ENV_VARS))
        if isinstance(refgenie_cfg_path, str) and not os.path.exists(refgenie_cfg_path):
            # file path not found, initialize a new config file
            print("File '{}' does not exist. Initializing refgenie genome configuration file.".format(refgenie_cfg_path))
            rgc = RGC(entries={CFG_FOLDER_KEY: os.path.dirname(refgenie_cfg_path)})
            rgc.initialize_config_file(filepath=refgenie_cfg_path)
        else:
            print("Reading refgenie genome configuration file from file: {}".format(refgenie_cfg_path))
            rgc = RGC(filepath=refgenie_cfg_path)
        try:
            # get path to the chrom.sizes asset
            chrom_sizes = rgc.get_asset(genome_name=args.genome, asset_name="fasta", tag_name="default",
                                        seek_key="chrom_sizes")
        except RefgenconfError:
            # if chrom.sizes not found, pull it first
            print("Could not determine path to chrom.sizes asset, pulling")
            rgc.pull_asset(genome=args.genome, asset="fasta", tag="default")
            chrom_sizes = rgc.get_asset(genome_name=args.genome, asset_name="fasta", tag_name="default",
                                        seek_key="chrom_sizes")
        print("Determined path to chrom.sizes asset: {}".format(chrom_sizes))
    else: 
        raise OSError("Could not determine path to chrom.sizes file. "
                    "Use --chrom-sizes argument to provide it, or "
                    "--rfg-config argument to  provide the path to genome config file and determine chrom.sizes asset via Refgenie.")
    
    return chrom_sizes


def main():
    pm = pypiper.PipelineManager(name="bedmaker", outfolder=logs_dir, args=args) # ArgParser and add_pypiper_args

    # Define target folder for converted files and implement the conversions; True=TF_Chipseq False=Histone_Chipseq

    print("Got input type: {}".format(args.input_type))
    if not args.input_type == "bed":
        print("Converting {} to BED format".format(args.input_file))

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
        input_file = os.path.join(os.path.dirname(args.output_bed),os.path.splitext(file_name)[0])
        with gzip.open(args.input_file, "rb") as f_in:
            with open(input_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        pm.clean_add(input_file)

    temp_bed_path = os.path.splitext(args.output_bed)[0]

    if args.input_type == "bedGraph":
        cmd = bedGraph_template.format(input=input_file, output=temp_bed_path, width=width)
    elif args.input_type == "bigWig":
        cmd = bigWig_template.format(input=input_file, output=temp_bed_path, width=width)
    elif args.input_type == "wig": 
        chrom_sizes = get_chrom_sizes()
        # define a target for temporary bw files
        temp_target = os.path.join(bed_parent, file_id + ".bw")
        pm.clean_add(temp_target)
        cmd1 = wig_template.format(input=input_file, intermediate_bw=temp_target, chrom_sizes=chrom_sizes, width=width)
        cmd2 = bigWig_template.format(input=temp_target, output=temp_bed_path, width=width)
        cmd = [cmd1, cmd2]
    elif args.input_type == "bigBed":
        cmd = bigBed_template.format(input=input_file, output=temp_bed_path)
    elif args.input_type == "bed":
        cmd = bed_template.format(input=args.input_file, output=args.output_bed)
    else:
        raise NotImplementedError("'{}' format is not supported".format(args.input_type))

    gzip_cmd = gzip_template.format(unzipped_converted_file=temp_bed_path)
    if args.input_type != "bed" or input_extension != ".gz":
        if not isinstance(cmd, list):
            cmd = [cmd]
        cmd.append(gzip_cmd)
    pm.run(cmd, target=args.output_bed)

    print("Generating bigBed files for {}".format(args.input_file))
    bedfile_name = os.path.split(args.output_bed)[1]
    fileid = os.path.splitext(os.path.splitext(bedfile_name)[0])[0]
    # Produce bigBed (bigNarrowPeak) file from peak file 
    bigNarrowPeak = os.path.join(args.output_bigbed, fileid + ".bigBed")
    if args.input_type != "bigBed":
        temp = tempfile.NamedTemporaryFile(dir=args.output_bigbed, delete=False)
        if not os.path.exists(bigNarrowPeak):            
            pm.clean_add(temp.name)
            cmd = ("zcat " + args.output_bed + "  | sort -k1,1 -k2,2n | awk '{print $1,$2,$3}' > " + temp.name)
            pm.run(cmd, temp.name, clean=True)

            chrom_sizes = get_chrom_sizes()

            cmd = ("bedToBigBed -type=bed3 " +
                    temp.name + " " + chrom_sizes + " " + bigNarrowPeak)
            pm.run(cmd, bigNarrowPeak, nofail=True)
    else:
        cmd = "ln -s {input} {output}".format(input=args.input_file, output=bigNarrowPeak)
        pm.run(cmd)
        
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
