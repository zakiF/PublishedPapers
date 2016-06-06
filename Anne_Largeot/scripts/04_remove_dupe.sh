# Filter, mark & remove duplicate in the merged bam files

#!/bin/sh

module load apps/samtools/0.1.19/gcc/4.4.7

source /data/acbb/to.keep/cwirth/troodonProcessing/troodonHeader.sh

working_dir=/data/stemcell/sczaki/project/anne/GL16-18

merged_bam=${working_dir}/bam_merged
merged_bam_list=(`ls ${merged_bam}/*.bam`)


logDir=${working_dir}/log/remove_dupe
mkdir ${logDir}

commandLog=${logDir}/dedduped_${startTime}.txt
touch ${commandLog}

jname=dedduped_${startTime}

# ----------------------------------------------------------------------------- #
# Remove Duplicates
# ----------------------------------------------------------------------------- #

temp_BAM_folder=${merged_bam}/dedduped_bam
mkdir ${temp_BAM_folder}

picardMetricsDir_temp=${temp_BAM_folder}/picard_metrics
mkdir ${picardMetricsDir_temp}

for file in ${merged_bam_list[*]}
do
        sampleName=`basename $file .bam`
        markoutput=${sampleName}_remove_dupe.bam

         waitForQueue

        command="java -jar /home/zfadlullah/prog/picard-tools-1.129/picard.jar MarkDuplicates INPUT=${file} OUTPUT=${temp_BAM_folder}/${markoutput} METRICS_FILE=${picardMetricsDir_temp}/${f}.txt REMOVE_DUPLICATES=TRUE"
        echo "${command}" >> ${commandLog}

  createQsubScript "${command}" ${referrer}
  /apps/modules/pkg/clusterware/torque/5.1.0-1/gcc-4.4.7/bin/qsub -N ${jname} -e ${logDir}/dedduped_${sampleName}_${startTime}_err.txt -o ${logDir}/dedduped_${sampleName}_${startTime}_out.txt ${referrer}

done
waitForName ${jname}

source /data/acbb/to.keep/cwirth/troodonProcessing/troodonFooter.sh