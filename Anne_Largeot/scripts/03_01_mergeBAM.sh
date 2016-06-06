#!/bin/sh
#
# Set name of job
#PBS -N mergeBAM

# Set the number of nodes and processors per node
#PBS -l nodes=2:ppn=8

# Define job queue: bioinfo (short), genome (long), vlm (high mem)
#PBS -q bioinf.q

# Mail alert at (b)eginning, (e)end and (a)bortion of execution
#PBS -m bea
#PBS -M zaki.fadlullah@cruk.manchester.ac.uk

# or you can define where to write output
#PBS -o /lustre/data/stemcell/sczaki/project/anne/GL16-18/log/merge
#PBS -e /lustre/data/stemcell/sczaki/project/anne/GL16-18/log/merge
#PBS -j oe

# Export all my environment variables to the job
#PBS -V

module load apps/R/3.1.0/gcc-4.4.7
module load apps/samtools/0.1.19/gcc/4.4.7


Rscript /home/zfadlullah/script/zaki/anne/GL16-18/03_02_mergeBAM.r






