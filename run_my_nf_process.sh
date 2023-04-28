#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cs
#SBATCH --time=5:00:00
#SBATCH --mem=5GB
#SBATCH --job-name=test_nf
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rj931@nyu.edu
#SBATCH --output=slurm_%j.out
#SBATCH --array=1-3


module purge

# loading the nextflow module below

module load nextflow/21.10.6

# giving the path to the txt file that contains the sample id, forward and reverse reads

PE_list='/scratch/work/courses/BI7653/hw2.2023/week2_fastqs.txt'

# using the array number to get the contents on that line, then separating the line by columns
# to get the name fo the foward read and name of the reverse read

line="$(less "$PE_list" | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)"
file_f="$(printf "%s" "${line}"| cut -f2)"
file_r="$(printf "%s" "${line}"| cut -f3)"

# below I realize i needed to put the name onto the location of the directory 
# this gets the full path of the foward file and reverse file

read_f='/scratch/work/courses/BI7653/hw2.2023/'$file_f
read_r='/scratch/work/courses/BI7653/hw2.2023/'$file_r

# later I need to change the name to from the input name to an output name used in fastp as the 
# output name later in the first process

name_f="$(basename $file_f .flit.fastq.gz).fastq.gz"
name_r="$(basename $file_r .flit.fastq.gz).fastq.gz"

# these echo lines are just for debugging

echo $line
echo $read_f
echo $read_r
echo $name_f
echo $name_r


# now call nextflow and tell it to run the name of my nextflow pipeline called final_project.nf
# I need to send all the variables i created above to the nf file so the pipeline can access
# them. i do this using the param parameters.
# I do not think i can use resume and could not figure out where to put it now that i have the 
# parameters there.

nextflow run final_project.nf --param1 $read_f --param2 $read_r --param3 $name_f \
--param4 $name_r
