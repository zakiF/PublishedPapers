
module load apps/samtools/0.1.19/gcc/4.4.7
module load apps/homer/4.7/gcc-4.4.7

working_dir=/lustre/data/stemcell/sczaki/project/anne

merged_bam=${working_dir}/bam_merged/chr_filtered/removed_dupe
merged_bam_list=(`ls ${merged_bam}/*.bam`)


for file in ${merged_bam_list[*]}
do
        f=`basename $file _remove_dupe.bam`
        makeTagDirectory ${merged_bam}/${f} ${file} -genome mm10 
done


#Proceed to make UCSC bigWig file


 f=`basename $file .bam`
 tagfiles=${merged_bam}${f}

for file in ${tagfiles}
do
	makeUCSCfile ${file} -o auto -noadj -bigWig /lustre/data/stemcell/sczaki/file/ref/mm10.chr.size -style chipseq -fsize 1e20
done

