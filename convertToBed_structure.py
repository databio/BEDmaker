#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os  
import sys  

parser = ArgumentParser(
	description = "A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("--input_file", help="path to the input file")
parser.add_argument("--Chip_exp", help="is it a ChiP-Seq TF experiment or a Histone modification ChiP-Seq experiment")
parser.add_argument("--input_type", help="a bigwig or a bedgraph file that will be converted into BED format") 
parser.add_argument("--outfolder", default="output", help=" folder to put the converted BED files in")

parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "common", "looper"], 
											args= ["input_file", "outfolder"],
											required=["input_file", "input_type", "Chip_exp"])

args = parser.parse_args() # all pipeline arguments should be defined by now 

if args.input_file is None:
	parser.print_help()
	sys.exit()

# get the file, exp, and file_type arguments for the .sh script 
file = args.input_file
exp = args.Chip_exp
file_type = args.input_type
outfolder = args.outfolder 

 # get the fileID argument for the .sh script 
 split_file_path = os.path.split(file) 

 file_tail =  split_file_path[1]
 dot_find = file_tail.find(".")
 if dot_find < 0 :
 	fileID = file_tail 
 else:
 	fileID = file_tail[0: dot_find]



# create Pipeline manager object" -------------------------------------------------------------------------------
 

pm = pypiper.PipelineManager(name="convertToBed_pipeline", outfolder=outfolder, args=args) #args defined with ArgParser and add_pypiper_args

command = "./file_conv_pseudocode.sh --input_file= %s --Chip_exp= %s --input_type= %s" % (file, exp, file_type)




