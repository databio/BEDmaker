#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os  

parser = ArgumentParser(description = "A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("-f", "--input-file", help="path to the input file", type=str)
parser.add_argument("-c", "--chip-exp", help="is it a ChIP-Seq TF experiment or a Histone modification ChiP-Seq experiment", type=str)
parser.add_argument("-t", "--input-type", help="a bigwig or a bedgraph file that will be converted into BED format")
parser.add_argument("-o", "--outfolder", default="output", help="folder to put the converted BED files in")

# parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "common", "looper"], args=["input_file", "outfolder"],
# 											required=["--input-file", "--input-type", "--chip-exp"])

args = parser.parse_args()  # all pipeline arguments should be defined by now

# commands templates
# bedGraph to bed
bed_graph_template = "macs2 bdgpeakcall -i {input} -o {output}"
# bigWig to bed
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


if args.input_type == "bedgraph":
    print("Converting {} to BED format".format(args.input_file))
    target = get_bed_path(args.input_file, args.outfolder)
    cmd = bed_graph_template.format(input=args.input_file, output=target)
else:
    raise NotImplementedError("Other conversions are not supported yet")

pm.run(cmd, target=target)
pm.stop_pipeline()
