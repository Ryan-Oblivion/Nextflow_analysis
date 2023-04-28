// define module


FASTP='fastp/intel/0.20.1'

BWA='bwa/intel/0.7.17'

SAMTOOLS='samtools/intel/1.14'

PICARD='picard/2.23.8'

// use picard or gatk?

GATK='gatk/4.3.0.0'

// this pipeline is created using the dsl2 commands
// so I can access that by using the next line
nextflow.enable.dsl=2



// first process below calles process_fqs_fastp

process process_fqs_fastp {



input:
// for the input i am using the parameters passed from the array .sh job but the names changed 
// here because of syntax seen in the workflow section at the end of this file

path read_f
path read_r
path out_n_f
path out_n_r


output:
// i used the name of the output created in the .sh slurm job and gave it a new variable 
// in the workflow section. then gave it a new variable for the next process to 
// call it with

path out_n_f, emit: fastp_file1
path out_n_r, emit: fastp_file2

// below is using unix coding to call the fastp module and give the inputs and outputs

"""
#!/bin/env bash
# if this works i can just do the slurm array here and get it to run 3 fastp arrays



module load $FASTP




fastp \
-i $read_f \
-I $read_r \
-o $out_n_f \
-O $out_n_r \
--length_required 76 \
--n_base_limit 50 \
--detect_adapter_for_pe \





"""

}

// this is the second process called align_bwa. the inputs are the files generated from the last 
// process and given the output channel names of fastp_file1 fastp_file2
// which can now be used as the input names


nextflow.enable.dsl=2


process align_bwa {

input: 
path fastp_file1
path fastp_file2
path index_file
path ref

output:
// this output creates a bam and gives it the bam_file name to channel it to the next process

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
// this has a new input that changes the name to a the samples name, so we have different 
// bam files at the end of each job. this name was made in the .sh file and was changed 
// to a new variable in the workflow

path bam_file
path out_n_f

output:
// for the output i added .bam to that sample name i wanted to be specific to each slurm array 
// job

path "${out_n_f}.bam", emit: coor_sorted_bam

"""


module load $PICARD
module load $SAMTOOLS

java -Xmx44g -jar \$PICARD_JAR SortSam \
I=$bam_file \
O="${out_n_f}.bam" \
SORT_ORDER=coordinate \

# we always have to create a BAM index file on any coordinate sorted BAM. 
# NOT possible to do so if it is not coordinate sorted

samtools index -b "${out_n_f}.bam" 

"""
}


// here is the workflow section that manages all the processes and channels
// the first four lines show how to get the inputs from the .sh file and give them new 
// variables that nextflow can keep track of.
// then you also get the index file and ref to have their own unique variables

// now we get to the main section where we tell nextflow wich process runs first second and
// third just by the order we place them in. then in the parameters of each process i 
// place the inputs i want that process to take. this is how nextflow manages its workflow

workflow {

read_f = Channel.fromPath(params.param1)
read_r = Channel.fromPath(params.param2)
out_n_f = Channel.fromPath(params.param3)
out_n_r = Channel.fromPath(params.param4)


index_file = Channel.fromPath(params.dir_ref_files)
ref = Channel.fromPath(params.ref)

main:
process_fqs_fastp(read_f, read_r, out_n_f, out_n_r)

align_bwa(process_fqs_fastp.out.fastp_file1, process_fqs_fastp.out.fastp_file2, index_file, ref)
coor_sort_picard(align_bwa.out.bam_file, out_n_f)


//coor_sort_picard.out.coor_sorted_index_bam.view()
//process_fqs_fastp.out.fastp_file1.view()
//process_fqs_fastp.out.fastp_file2.view()
//align_bwa.out.bam_file.view()

}




