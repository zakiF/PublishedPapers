
# Extend the and merge the range in bed files

# 10/09/2015

module load apps/bedtools/2.20.1/gcc-4.4.7 


workDir=/data/stemcell/sczaki/project/anne/MT2_MOZ/MT2/p_value_0.01/clean/UCSC_edited
bed_list=(`ls ${workDir}/*.bed`)

# Change these parameters
# 1) The range (line 16)
# 2) Name of extension (line 28)


# The number of bases in bp to extend the bed
range=250

# mm10 genome chr size
knownSize=/scratch/zaki/file/ref/mm10.chr.size_knownChr

# Create a folder to store the output
extendDir=${workDir}/Extended_bed
mkdir ${extendDir}

for file in ${bed_list[*]}
do
        f=`basename $file _peaks_narrowPeak_clean.bed_UCSC.bed` # Change this line
        edited=${f}_extended_${range}.bed
        bedtools slop -i ${file} -b ${range} -g ${knownSize} | bedtools merge -i - > ${extendDir}/${edited}
done



