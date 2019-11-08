#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os  

parser = ArgumentParser(description = "A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("-f", "--input-file", help="path to the input file", type=str)
parser.add_argument("-c", "--chip-exp", help="is it a ChIP-Seq TF experiment or a Histone modification ChiP-Seq experiment", type=str)
parser.add_argument("-t", "--input-type", help="a bigwig or a bedgraph file that will be converted into BED format")
parser.add_argument("-o", "--outfolder", default="output", help="folder to put the converted BED files in")


# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                            required=["--input-file", "--input-type", "--chip-exp", "--outfolder"])

args = parser.parse_args()  

# commands templates
# bedGraph to bed
bedGraph_template = "macs2 bdgpeakcall -i {input} -o {output}"
# bigWig to bed
bigwig_template = "bigWigToBedgraph -i {input} | macs2 bdgpeakcall -i {intermediate}  -o {output}"
# big_wig_cmd = ""

pm = pypiper.PipelineManager(name="bed_maker", outfolder=args.outfolder, args=args) #args defined with ArgParser and add_pypiper_args


def get_bed_path(current_path, outfolder):
    """
    Swap the file extension and change the directory

    :param str current_path: current path to the file to be converted
    :param str outfolder: output directory to place the file with the swapped extension in
    :return str: path to the file with swapped extension
    """
    file_name = os.path.basename(current_path)
    file_id = os.path.splitext(file_name)[0]
    return os.path.join(outfolder, file_id + ".bed")

# Define target folder for converted files
target = get_bed_path(args.input_file, args.outfolder)


if args.chip_exp == True:
    if args.input_type == "bedGraph":
        print("Converting {} to BED format".format(args.input_file))
        cmd = bedGraph_template.format(input=args.input_file, output=target)
    elif args.input_type == "bigWig":
        print("Converting {} to BED format".format(args.input_file))
        cmd = bigwig_template.format(input=args.input_file, 
                                    intermediate=os.path.splitext(os.path.basename(args.input_file))[0] + ".bedgraph",  
                                    output=target)
else:
    raise NotImplementedError("Other conversions are not supported yet")

pm.run(cmd, target=target)
pm.stop_pipeline()