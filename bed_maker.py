#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os
import sys
import pyBigWig
from refgenconf import RefGenConf as RGC, select_genome_config, RefgenconfError, CFG_ENV_VARS, CFG_FOLDER_KEY
import yacman

parser = ArgumentParser(description="A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("-f", "--input-file", help="path to the input file", type=str)
parser.add_argument("-c", "--chip-exp", help="is it a ChIP-Seq TF experiment or a Histone modification ChiP-Seq experiment", type=bool)
parser.add_argument("-t", "--input-type", help="a bigwig or a bedgraph file that will be converted into BED format")
parser.add_argument("-g", "--genome", help="reference genome")
parser.add_argument("-r", "--rfg-config", help="file path to the genome config file", type=str)
parser.add_argument("-o", "--output-file", help="path to the output BED files", type=str)
#parser.add_argument("-s", "--sample-name", help="name of the sample used to systematically build the output name", type=str)


# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                  required=["--input-file", "--input-type", "--chip-exp"])

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


# SET OUTPUT FOLDERS
# use output parent argument from looper 
out_parent = args.output_parent

file_name = os.path.basename(args.input_file)
file_id = os.path.splitext(file_name)[0]
sample_folder = os.path.join(out_parent, file_id) #specific output folder for each sample log and stats

bed_parent = os.path.dirname(args.output_file)
if not os.path.exists(bed_parent):
    print("Output directory does not exist. Creating: {}".format(bed_parent))
    os.makedirs(bed_parent)


def main():
    pm = pypiper.PipelineManager(name="bed_maker", outfolder=sample_folder, args=args) # ArgParser and add_pypiper_args

    # Define target folder for converted files and implement the conversions; True=TF_Chipseq False=Histone_Chipseq

    print("Got input type: {}".format(args.input_type))
    print("Converting {} to BED format".format(args.input_file))

    # Define whether Chip-seq data has broad or narrow peaks
    width = "bdgbroadcall" if not args.chip_exp else "bdgpeakcall"

    # Call pyBigWig to ensure bigWig and bigBed files have the correct format
    if args.input_type in ["bigWig", "bigBed"]:
        big_check = pyBigWig.open(args.input_file)

    # if args.chip_exp:
    if args.input_type == "bedGraph":
        cmd = bedGraph_template.format(input=args.input_file, output=args.output_file, width=width)
    elif args.input_type == "bigWig" and big_check.isBigWig():
        cmd = bigWig_template.format(input=args.input_file, output=args.output_file, width=width)
    elif args.input_type == "wig": 
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
        # define a target for temporary bw files
        temp_target = os.path.join(sample_folder, file_id + ".bw")
        pm.clean_add(temp_target)
        cmd1 = wig_template.format(input=args.input_file, intermediate_bw=temp_target, chrom_sizes=chrom_sizes, width=width)
        cmd2 = bigWig_template.format(input=temp_target, output=args.output_file, width=width)
        cmd = [cmd1, cmd2]
    elif args.input_type == "bigBed" and big_check.isBigBed():
        cmd = bigBed_template.format(input=args.input_file, output=args.output_file)
    elif args.input_type == "bed":
        cmd = bed_template.format(input=args.input_file, output=args.output_file)
    else:
        raise NotImplementedError("'{}' format is not supported".format(args.input_type))

    pm.run(cmd, target=args.output_file)
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
