#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os
import sys
import pyBigWig

parser = ArgumentParser(description = "A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("-f", "--input-file", help="path to the input file", type=str)
parser.add_argument("-c", "--chip-exp", help="is it a ChIP-Seq TF experiment or a Histone modification ChiP-Seq experiment", type=bool)
parser.add_argument("-t", "--input-type", help="a bigwig or a bedgraph file that will be converted into BED format")
parser.add_argument("-g", "--genome", help="reference genome")
#parser.add_argument("-o", "--outfolder", default="output", help="folder to put the converted BED files in")


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
wig_template =  "wigToBigWig {input} {chrom_sizes} /dev/stdout | bigWigToBedGraph /dev/stdin  /dev/stdout | macs2 {width} -i /dev/stdin -o {output}"

out_parent = args.output_parent # use output parent argument from looper 

# def get_bed_path(current_path, outfolder):
#     """
#     Swap the file extension and change the directory

#     :param str current_path: current path to the file to be converted
#     :param str outfolder: output directory to place the file with the swapped extension in
#     :return str: path to the file with swapped extension
#     """
#     file_name = os.path.basename(current_path)
#     file_id = os.path.splitext(file_name)[0]
#     return os.path.join(outfolder, file_id + ".bed")

file_name = os.path.basename(args.input_file)
file_id = os.path.splitext(file_name)[0]
sample_folder = os.path.join(out_parent, file_id)

def main():
    pm = pypiper.PipelineManager(name="bed_maker", outfolder=sample_folder, args=args) # ArgParser and add_pypiper_args

    # Define target folder for converted files and implement the conversions; True=TF_Chipseq False=Histone_Chipseq
    #target = get_bed_path(args.input_file, outfolder2)
    target = os.path.join(sample_folder, file_id + ".bed")

    print("Got input type: {}".format(args.input_type))
    print("Converting {} to BED format".format(args.input_file))

    # Define whether Chip-seq data has broad or narrow peaks
    width = "bdgbroadcall" if not args.chip_exp else "bdgpeakcall"
    
    
    # Call pyBigWig to ensure bigWig and bigBed files have the correct format
    if args.input_type == "bigWig" or args.input_type == "bigBed":
        big_check = pyBigWig.open(args.input_file)
    
    
    # if args.chip_exp:
    if args.input_type == "bedGraph":
        cmd = bedGraph_template.format(input=args.input_file, output=target, width=width)
    elif args.input_type == "bigWig" and big_check.isBigWig() :
        cmd = bigWig_template.format(input=args.input_file, output=target, width=width)
    elif args.input_type == "wig": 
        if args.genome == "hg19":
            cmd = wig_template.format(input=args.input_file, output=target, chrom_sizes="ftp://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/hg19.chrom.sizes", width=width)
        elif args.genome == "hg38":
             cmd = wig_template.format(input=args.input_file, output=target, chrom_sizes=hg38.chrom.sizes, width=width)
    elif args.input_type == "bigBed" and big_check.isBigBed():
        cmd = bigBed_template.format(input=args.input_file, output=target)
    else:
        raise NotImplementedError("Other conversions are not supported")

    pm.run(cmd, target=target)
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
