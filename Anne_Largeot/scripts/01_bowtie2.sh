### Align GL16 (Anne-ChipSeq) sequencing data using bowtie2 to mouse genome assembly mm10 using bowtie2 (version 2.2.1)

module load apps/bowtie2/2.2.1/gcc/4.4.7
module load apps/samtools/0.1.19/gcc/4.4.7 

source /lustre/data/acbb/to.keep/cwirth/troodonProcessing/troodonHeader.sh

#Location of bowtie2inex
#Using mm10 UCSC reference
bowtie2Index=/lustre/scratch/zaki/file/ref/bowtie2/mm10/mm10
#Location of seuqencing files by lane
fastqDir=/lustre/facilities/sequencers/seq.02_fromHPC1/140224_D00128_0110_BH8FAAADXX_Unaligned/Unaligned/Project_DefaultProject
r1files=(`ls ${fastqDir}/*/GL16*.gz`)

#Creating a directory to store all output
projectOutputDir=/lustre/data/stemcell/sczaki/project/anne/GL16-18
mkdir ${projectOutputDir}	
	
#Directory for samfiles
samDir=${projectOutputDir}/sam
mkdir ${samDir}

#Directory for output logs
logDir=${projectOutputDir}/log
mkdir ${logDir}

#Mapping stats
mappingStatsDir=${projectOutputDir}/bowtie2_mapping_statistics
mkdir ${mappingStatsDir}

#Recording command used for bowtie2mapping
commandLog=${logDir}/bowtie2Commands_${startTime}.txt
touch ${commandLog}

jname=bowtie2_${startTime}


#Creating a loop function to conduct allignment on all paired fastq file in the working directory described
for inputFastq1 in ${r1files[*]}; do
# %% is a pattern matching command - http://tldp.org/LDP/abs/html/parameter-substitution.html
  noExt=${inputFastq1%%_R1_001.fastq.gz}
  
  sampleName=$(basename ${noExt})
  output=${samDir}/${sampleName}.sam

  waitForQueue


  command="bowtie2 -p 1 -q -x ${bowtie2Index} -U ${inputFastq1}  -S ${output}"
  echo "${command}" >> ${commandLog}


  createQsubScript "${command}" ${referrer}
  qsub -N ${jname} -e ${mappingStatsDir}/${sampleName}_mappingStats.txt -o ${logDir}/bowtie2_${sampleName}_${startTime}_out.txt ${referrer}

done
waitForName ${jname}

###################################What does this script do
source /lustre/data/acbb/to.keep/cwirth/troodonProcessing/troodonFooter.sh

