// define module


FASTP='fastp/intel/0.20.1'

BWA='bwa/intel/0.7.17'

SAMTOOLS='samtools/intel/1.14'

PICARD='picard/2.23.8'

// use picard or gatk?

GATK='gatk/4.3.0.0'

nextflow.enable.dsl=2
//


params.input1 = './SRR708363_1.filt.fastq.gz'
params.input2 = './SRR708363_2.filt.fastq.gz'

process process_fqs_fastp {

input:
path reads1
path reads2

output:
path 'SRR708363_1.fastq.gz', emit: fastp_file1
path 'SRR708363_2.fastq.gz', emit: fastp_file2
 


"""
#!/bin/env bash
# if this works i can just do the slurm array here and get it to run 3 fastp arrays



module load $FASTP





fastp \
-i $reads1 \
-I $reads2 \
-o SRR708363_1.fastq.gz \
-O SRR708363_2.fastq.gz \
--length_required 76 \
--n_base_limit 50 \
--detect_adapter_for_pe \





"""

}

// this is the second process called align_bwa. the inputs are the files generated from the last 
// process and given the output channel names of fastp_file1 fastp_file2


nextflow.enable.dsl=2


process align_bwa {

input: 
path fastp_file1
path fastp_file2
path index_file
path ref

output:
path 'aligned_reads_w_header.bam', emit: bam_file


"""
#!/bin/env bash

module load $BWA
module load $SAMTOOLS


bwa index -a bwtsw $ref

bwa mem -M -t 8 $ref $fastp_file1 $fastp_file2 > aligned_reads.sam

samtools view -b -h aligned_reads.sam -o aligned_reads_w_header.bam

"""


}

// new process to coordinate sort the bam file

process coor_sort_picard{
cpus 10
memory '46 GB'
executor 'slurm'

input:
path bam_file

output:
path 'First_file_coor_sorted.bam', emit: coor_sorted_bam

"""


module load $PICARD
module load $SAMTOOLS

java -Xmx44g -jar \$PICARD_JAR SortSam \
I=$bam_file \
O=First_file_coor_sorted.bam \
SORT_ORDER=coordinate \

# we always have to create a BAM index file on any coordinate sorted BAM. 
# NOT possible to do so if it is not coordinate sorted

samtools index -b First_file_coor_sorted.bam 

"""
}


workflow {

reads1 = Channel.fromPath(params.input1)
reads2 = Channel.fromPath(params.input2)
index_file = Channel.fromPath(params.dir_ref_files)
ref = Channel.fromPath(params.ref)

main:
process_fqs_fastp(reads1, reads2)
align_bwa(process_fqs_fastp.out.fastp_file1, process_fqs_fastp.out.fastp_file2, index_file, ref)
coor_sort_picard(align_bwa.out.bam_file)


//coor_sort_picard.out.coor_sorted_index_bam.view()
//process_fqs_fastp.out.fastp_file1.view()
//process_fqs_fastp.out.fastp_file2.view()
//align_bwa.out.bam_file.view()

}




