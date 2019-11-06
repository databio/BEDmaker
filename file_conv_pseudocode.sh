#!/usr/bin/bash

if [exp == "TF"]
	then
		if [ file_type == "bedGraph" ] 
			macs2 bdgpeakcall -t file 
		elif [ file_type == "bigWig" ]
			bigWigToBedgraph file "input_file.bedgraph" | macs2 bdgpeakcall -t input_file.bedgraph
elif [ exp == "Histone_mark" ]
	then
		if [ file_type == "bedgraph"]

			






