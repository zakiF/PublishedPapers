
# Combine any two peaks which are close in proximity ~100bp

# 18/09/2015

module load apps/bedtools/2.20.1/gcc-4.4.7 


workDir=/data/stemcell/sczaki/project/anne/MT2_MOZ/MT2/p_value_0.01/clean/UCSC_edited/Extended_bed
bed_list=(`ls ${workDir}/*.bed`)

# Change these parameters
# 1) The range (line 16)
# 2) Name of extension (line 28)


# The number of bases in bp to combine if two peaks are close
range=100

# mm10 genome chr size
knownSize=/scratch/zaki/file/ref/mm10.chr.size_knownChr

# Create a folder to store the output
mergeDir=${workDir}/merge
mkdir ${mergeDir}

for file in ${bed_list[*]}
do
        f=`basename $file .bed` # Change this line
        edited=${f}_M.bed
        bedtools merge -d ${range} -i ${file} > ${mergeDir}/${edited}
done


# ----------------------------------------------------------------------------- #
# Generate a bigBed for these bed
# ----------------------------------------------------------------------------- #

# Location of bed2bigBed program

bedToBigBed=/home/zfadlullah/bin/bedToBigBed

bed_list=(`ls ${mergeDir}/*.bed`)
for file in ${bed_list[*]}
do
        f=`basename $file .bed`
        edited=${f}.bb
        ${bedToBigBed} ${file} ${knownSize} ${mergeDir}/${edited}
done



