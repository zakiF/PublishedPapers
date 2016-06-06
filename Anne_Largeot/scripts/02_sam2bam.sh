### Convert sam to bam (lane-wise)

module load apps/samtools/0.1.19/gcc/4.4.7 

source /lustre/data/acbb/to.keep/cwirth/troodonProcessing/troodonHeader.sh

#Make sure its the same output directory
projectOutputDir=/lustre/data/stemcell/sczaki/project/anne/GL16-18

#Output of mapping with bowtie2 should be in this folder
samDir=${projectOutputDir}/sam
infiles=(`ls ${samDir}/*.sam`)

#Desired location of bam file
bamDir=${projectOutputDir}/bam_lanewise
mkdir ${bamDir}

#Desired location for temporary bam file (to calculate lane wise duplication level)
bamDir2=${projectOutputDir}/bam_lanewise_dupe_temp
mkdir ${bamDir2}

logDir=${projectOutputDir}/log

commandLog=${logDir}/sam2bamCommands_${startTime}.txt
touch ${commandLog}

jname=sam2bam_${startTime}

for samInput in ${infiles[*]}; do

  sampleName=$(basename ${samInput} .sam)
  output=${bamDir}/${sampleName}.bam
  output2=${bamDir2}/${sampleName}.bam

  waitForQueue

#This is the BAM file to be used for downstream purpose
command="samtools view -bS -q 10 -F 260 ${samInput} -o ${output}"
# The -q 10 skip alignments with MAPQ score less than 10 (most likely to be multimapping reads)
# The -F 260 flag, removes reads that are :
# Unmapped
# Not primary alignment
  echo "${command}" >> ${commandLog}

  createQsubScript "${command}" ${referrer}
  qsub -N ${jname} -e ${logDir}/sam2bam_${sampleName}_${startTime}_err.txt -o ${logDir}/sam2bam_${sampleName}_${startTime}_out.txt ${referrer}

#For mark duplicates and getting coverage calculation
command_2="samtools view -bS -F 260 ${samInput} -o ${output2}"
# The -F 260 flag, removes reads that are :
# Unmapped
# Not primary alignment
  #No need to record commands, as the bam will only act as temporary files
  #echo "${command_2}" >> ${commandLog}

  createQsubScript "${command_2}" ${referrer}
  qsub -N ${jname} -e ${logDir}/sam2bam_dupe_temp_${sampleName}_${startTime}_err.txt -o ${logDir}/sam2bam_dupe_temp_${sampleName}_${startTime}_out.txt ${referrer}


done
waitForName ${jname}

source /lustre/data/acbb/to.keep/cwirth/troodonProcessing/troodonFooter.sh

