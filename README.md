# bed_maker
Pypiper pipeline to convert bigWig, bedGraph, bigBed and other types of files into BED format.

### Before running the pipeline: 

Users will need to format metadata as a Portable Encapsulated Project (PEP). Looper will use the PEP annotation sheet to run the pipeline on all the samples.  
  

### To run the pipeline: 

  ` looper run example/cfg.yaml`
 
### Dependencies:

pyBigWig

macs2

bigWigToBedGraph
