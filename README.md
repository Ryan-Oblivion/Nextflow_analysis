# Nextflow_analysis
Here I create a nextflow pipeline that takes pair end reads and outputs a coordinate sorted bam and bai file. Look at instructions file

The PE_read_pipeline.nf file contains 3 processes that take pair end (PE) reads and outputs, at the very end a coordinate sorted bam file and an index file (bai). Tools used are fastp, bwa mem, samtools, and picard.

The nextflow.config file contains some variables that will be used in the pipeline and it tells nextflow what executor to use. In this case it would be slurm. The config file also specifies how many resources to allocate, however in the last process I needed to specifically specify there how many resources.
I will update the pipeline so it can take many other inputs also without having to specify and hard code the primary input names into the PE_read_pipeline.nf 
