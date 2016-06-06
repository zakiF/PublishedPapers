#MACS2 modify bed


# Input files is bed file from MACS2

# ----------------------------------------------------------------------------- #
# 1) Remove black list position from .broadPeak file
# ----------------------------------------------------------------------------- #
#	Chr Start End maxScore
# Modification for annotation pupose in R package ChIPseeker

module load apps/bedtools/2.20.1/gcc-4.4.7

working_dir=/lustre/data/stemcell/sczaki/project/anne
bed_folder=${working_dir}/bam_merged/chr_filtered/removed_dupe/MACS2
bed_list=(`ls ${bed_folder}/*.broadPeak`)


oriClean_folder=${bed_folder}/clean
mkdir ${oriClean_folder}

# mm10 blacklisted region obtained using liftOver tools from mm9
# mm9 black list position obtained from 
# http://www.broadinstitute.org/~anshul/projects/mouse/blacklist/
# The following command was used to liftOever mm9 -> mm10
# ./liftOver /lustre/scratch/zaki/file/ref/blacklist/mm9-blacklist.bed /lustre/scratch/zaki/file/ref/blacklist/mm9ToMm10.over.chain.gz mm10-blacklist.bed unmapped_mm10-blacklist.bed


blacklist=/lustre/scratch/zaki/file/ref/blacklist/mm10-blacklist.bed

for file in ${bed_list[*]}
do
	f=`basename ${file} .broadPeak`
	edited=${f}_broadPeak_clean.bed
	bedtools intersect -a ${file} -b ${blacklist} -v > ${oriClean_folder}/${edited}
done


# ----------------------------------------------------------------------------- #
# 2) Remove blacklisted position in oringal xls output as well
# ----------------------------------------------------------------------------- #

xls_list=(`ls ${bed_folder}/*.xls`)

# Remove "#" symbols from .xls files
for file in ${xls_list[*]}
do
	f=`basename ${file} .xls`
	edited=${f}_edited.xls
	cat ${file} | awk '($0 !~ /#/)' |  tail -n +3 > ${oriClean_folder}/${edited}
done

# Filter out blacklisted position
ori_list=(`ls ${oriClean_folder}/*.xls`)
blacklist=/lustre/scratch/zaki/file/ref/blacklist/mm10-blacklist.bed

for file in ${ori_list[*]}
do
	f=`basename ${file} .xls`
	clean=${f}_clean.xls
	bedtools intersect -a ${file} -b ${blacklist} -v > ${oriClean_folder}/${clean}
done

# Remove the 1st phase edited bed file
rm ${oriClean_folder}/*_edited.xls




