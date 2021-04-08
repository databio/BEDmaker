#!/usr/bin/env python

from argparse import ArgumentParser
import pypiper
import os
import sys
import tempfile
import pandas as pd
import gzip
import shutil
from refgenconf import RefGenConf as RGC, select_genome_config, RefgenconfError, CFG_ENV_VARS, CFG_FOLDER_KEY
from yacman.exceptions import UndefinedAliasError

parser = ArgumentParser(description="A pipeline to convert bigwig or bedgraph files into bed format")


parser.add_argument("-f", "--input-file", help="path to the input file", type=str)
parser.add_argument("-n", "--narrowpeak", help="whether the regions are narrow (transcription factor implies narrow, histone mark implies broad peaks)", type=bool)
parser.add_argument("-t", "--input-type", help="a bigwig or a bedgraph file that will be converted into BED format", type=str)
parser.add_argument("-g", "--genome", help="reference genome", type=str)
parser.add_argument("-r", "--rfg-config", help="file path to the genome config file", type=str)
parser.add_argument("-o", "--output-bed", help="path to the output BED files", type=str)
parser.add_argument("--output-bigbed", help="path to the output bigBed files", type=str)
parser.add_argument("-s", "--sample-name", help="name of the sample used to systematically build the output name", type=str)
#parser.add_argument("--chrom-sizes", help="a full path to the chrom.sizes required for the bedtobigbed conversion", type=str)

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
if args.input_type != "bed":
    if input_extension == ".gz":
        output_bed = os.path.splitext(os.path.splitext(args.output_bed)[0])[0] + '.bed.gz'
    else:
        output_bed = os.path.splitext(args.output_bed)[0] + '.bed.gz'
else:
    if input_extension != ".gz":
        output_bed = args.output_bed + '.gz'
    else:
        output_bed = args.output_bed

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

pm = pypiper.PipelineManager(name="bedmaker", outfolder=logs_dir, args=args) # ArgParser and add_pypiper_args

def get_chrom_sizes():
    """
    Get chrom.sizes file path by input arg, or Refegenie.

    :return str: chrom.sizes file path
    """

    print("Determining path to chrom.sizes asset via Refgenie.")
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
        chrom_sizes = rgc.seek(genome_name=args.genome, asset_name="fasta", tag_name="default", seek_key='chrom_sizes')
        print (chrom_sizes)
    except (UndefinedAliasError, RefgenconfError):
        # if chrom.sizes not found, pull it first
        print("Could not determine path to chrom.sizes asset, pulling")
        rgc.pull(genome=args.genome, asset="fasta", tag="default")
        chrom_sizes = rgc.seek(genome_name=args.genome, asset_name="fasta", tag_name="default", seek_key='chrom_sizes')

    print("Determined path to chrom.sizes asset: {}".format(chrom_sizes))
    
    return chrom_sizes


def validate_genome_assembly(chrom_sizes, bed):
    """
    validate the regions in the input bed file matches the chrom.sizes

    :return bool
    """
    
    print("Validating chromsome numbers for {}.".format(output_bed))
    temp = os.path.join(args.output_bigbed, next(tempfile._get_candidate_names()))
    pm.clean_add(temp)

    cmd = ("cut -d' ' -f1 " + bed + " | sort -u > " + temp)
    pm.run(cmd, temp)
    cmd = ("grep -owf " + temp + " " + chrom_sizes + " | uniq | grep -vf - " + temp)
    chrom_num = pm.run(cmd, lock_name = "lock."+next(tempfile._get_candidate_names()))
    print (chrom_num)


    # df_cs = pd.read_csv(chrom_sizes, sep="\t", header=None)
    # df_bed = pd.read_csv(bed, sep=" ", header=None)

    # print("Validating chromsome numbers for {}.".format(output_bed))
    # if df_bed[0].isin(df_cs[0]).all(axis=None):
    #     print("Validating region coordinates for {}.".format(output_bed))
    #     out_of_range = pd.DataFrame()
    #     for index, row in df_bed.iterrows():
    #         size = df_cs[df_cs[0]==row[0]][1].values[0]
    #         if row[1] > size or row[2] > size:
    #             out_of_range = out_of_range.append(df_bed[index])
    #     if out_of_range.empty:
    #         return True
    #     else:
    #         print (out_of_range)
    #         print("The following regions in {} is out of range: {}".format(output_bed, out_of_range))
    #         return False

    # else:
    #     chrom_list = list(set(df_bed[0]).difference(df_cs[0]))
    #     print("{} are not found in the chrom.sizes file.".format(chrom_list))
    #     return False


def main():
    # pm = pypiper.PipelineManager(name="bedmaker", outfolder=logs_dir, args=args) # ArgParser and add_pypiper_args

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
        input_file = os.path.join(os.path.dirname(args.input_file),os.path.splitext(file_name)[0])
        with gzip.open(args.input_file, "rb") as f_in:
            with open(input_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        pm.clean_add(input_file)
        
    temp_bed_path = os.path.splitext(output_bed)[0]

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
        cmd = bed_template.format(input=input_file, output=temp_bed_path)
    else:
        raise NotImplementedError("'{}' format is not supported".format(args.input_type))

    if args.input_type == "bed" and input_extension != ".gz":
        gzip_cmd = gzip_template.format(unzipped_converted_file=input_file)
    else:
        gzip_cmd = gzip_template.format(unzipped_converted_file=temp_bed_path)

    if args.input_type != "bed" or input_extension != ".gz":
        if not isinstance(cmd, list):
            cmd = [cmd]
        if args.input_type == "bed":
            cmd.insert(0, gzip_cmd)
        else:
            cmd.append(gzip_cmd)
    pm.run(cmd, target=output_bed)

    print("Generating bigBed files for {}".format(args.input_file))
    
    bedfile_name = os.path.split(output_bed)[1]
    fileid = os.path.splitext(os.path.splitext(bedfile_name)[0])[0]
    # Produce bigBed (bigNarrowPeak) file from peak file 
    bigNarrowPeak = os.path.join(args.output_bigbed, fileid + ".bigBed")
    chrom_sizes = get_chrom_sizes()
    temp = os.path.join(args.output_bigbed, next(tempfile._get_candidate_names())) 
    
    if not os.path.exists(bigNarrowPeak):            
        pm.clean_add(temp)
        cmd = ("zcat " + output_bed + "  | awk '{print $1,$2,$3}' |  sort -k1,1 -k2,2n > " + temp)
        pm.run(cmd, temp)
        try:
            cmd = ("bedToBigBed " +
                    temp + " " + chrom_sizes + " " + bigNarrowPeak)
            pm.run(cmd, bigNarrowPeak)
        except:
            print("Fail to generating bigBed files for {}".format(args.input_file))
    
    # if args.input_type != "bigBed":temp_bed_path = os.path.splitext(output_bed)[0]
    #     chrom_sizes = get_chrom_sizes()
    #     temp = os.path.join(args.output_bigbed, next(tempfile._get_candidate_names())) 
    #     if not os.path.exists(bigNarrowPeak):            
    #         pm.clean_add(temp)
    #         cmd = ("zcat " + output_bed + "  | awk '{print $1,$2,$3}' |  sort -k1,1 -k2,2n > " + temp)
    #         pm.run(cmd, temp)
    #         if validate_genome_assembly(chrom_sizes, temp):
    #             cmd = ("bedToBigBed " +
    #                     temp + " " + chrom_sizes + " " + bigNarrowPeak)
    #             pm.run(cmd, bigNarrowPeak, nofail=True)
    #         else: 
    #             print("Fail to generating bigBed files for {}".format(args.input_file))     
    # else:
    #     cmd = "ln -s {input} {output}".format(input=args.input_file, output=bigNarrowPeak)
    #     pm.run(cmd)
        
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
