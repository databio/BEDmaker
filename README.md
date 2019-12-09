# bed_maker
Pypiper pipeline to convert bigWig, bedGraph, bigBed and other types of files into BED format.

### Before running the pipeline: 

Users will need to format metadata as a Portable Encapsulated Project (PEP). Looper will use the PEP yaml config file to run the pipeline on all the dataset samples.  
  
### Current formats supported:
- bedGraph
- bigBed
- bigWig
- wig

### To run the pipeline: 

  ` looper run example/cfg.yaml`
 
