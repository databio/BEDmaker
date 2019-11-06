#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os  
import sys  

parser = ArgumentParser(
	description = "A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("input_file", help="path to the input file")
parser.add_argument("--Chip_exp", help="is it a ChiP-Seq TF experiment or a Histone modification ChiP-Seq experiment")
parser.add_argument("--input_type", help="a bigwig or a bedgraph file that will be converted into BED format") 
parser.add_argument("--outfolder", default="output", help=" folder to put the BED files in")

parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "common", "looper"], 
											args= ["input_file", "output_folder"],
											required=["input_file", "input_type", "Chip_experiment"])
args = parser.parse_args()

 



# Wait to create Pipeline manager object" -------------------------------------------------------------------------------
outfolder = os.path.join() #define this later

pm = pypiper.PipelineManager(name="convertToBed_pipeline", outfolder=outfolder, args=args) #args defined with ArgParser and add_pypiper_args

command = 




