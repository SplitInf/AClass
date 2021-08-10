# AClass

# 0.1.2-beta	Added support for Ho et al. colour scheme
		Added function to specify output path
		training_memberships_link renamed to training_memberships_path
		Changed file output naming scheme. No longer have numeric value before filenames.
		prefix is now optional
		Various bug fixes

# 0.1.1-beta	Added support for R 4.1.0 and limma 3.48.1. Various bug fixes
		Added batch.process.raw()
		Changed nano.load() behavior to not load recurrsively to prevent mixing of batches 

# 0.1.0-beta	Encompass the following changes from initial release	
# version 200117  Added feature to allow SampleContent option in process.raw() and nano.norm().
#                 Adding feature to allow include_samples_path and omit_samples_path options in classify.data() 
#                 Fixed nano.norm(), nano.prep() and nano.plot() to handle single sample.
#                 Fixed nano.MDS() to exit gracefully when data has less than 3 samples.
#                 Fixed "3_",prefix,"_testing_summary_full.txt" format
#                 Changed condition in nano.prep to allow zero variance sample and to proceed and issue warning.
#                 Changed nano.load() to allow files with the keyword “.*NormalizedData.*.csv” to be processed
#
# version 190709  First release


#Usage

keep_file_path - *optional* files containing sample names (one sample per row, no header) to include for process.raw() and classify.data(). If this option is leave blank, all samples in raw_path will be loaded. Works in opposite way as omit_file_path
omit_file_path - *optional*
probes_rank_path - *optional* probe ranking list used for training. By default uses "probes_list.txt" 



#FAQ
Q: How should I process multiple batches / Concerns for batch effect.
A: It is recommended to process each NanoString batch separetly to minimize the effect of batch effect. It is known that poor quality samples tend to skew the overall normalization. nano.prep() would issue a warning when high variance outlier is detected. A run with majority failed samples would severely affect outcome and it is recomended to re-run without the failed samples.


Q: How to handle single sample analysis?
A: Raw data from the training set can be used to normalize with the single sample. 

"Training_data" directory contains:
Trainingdata.list - list of sample names used in pre-trained models
Training_NormalizedData.csv - Raw training data used in pre-trained models

The training samples to be placed together in the same folder as the single sample with the following workflow:

1.  Include the pre-trained models sample Training_NormalizedData.csv in the folder that contain the single new sample.
2.  Run process.raw() normally
3.  Include omit_file_path to "Trainingdata.list" in classify.data() to filter out the pre-trained models samples from the output
