# Merge and index BAM files in GL16-18
#
# 28 March 2015, Zaki


#working.dir <- "/lustre/data/compbio/hsleong/Alignments/GL22"
working.dir <- "/data/stemcell/sczaki/project/anne/GL16-18"
setwd(working.dir)

#dest.dir <- "/lustre/data/compbio/hsleong/Alignments/GL22/bam_merged"
dest.dir <- paste0(working.dir,"/bam_merged")
dir.create(dest.dir)

#input.dir <- "/lustre/data/compbio/hsleong/Alignments/GL22/bam_lanewise"
input.dir <- paste0(working.dir,"/bam_lanewise")
myfiles <- list.files(input.dir, pattern=".bam$")


#commands.dir <- "/lustre/data/compbio/hsleong/Alignments/GL22/log"
commands.dir <- paste0(working.dir,"/log/merge")
dir.create(commands.dir)
command.file <- file.path(commands.dir, "mergeBamCommands.txt")

  
output.dir <- "alignmentStats" # Create folder to store all results
dir.create(output.dir)

pnames <- gsub("_L00(1|2).bam", "", myfiles)
treatment <- vector("character", length(pnames))
treatment[grep("-Input-WT_", pnames)] <- "GL16_Input-WT"
treatment[grep("-Input-WT_", pnames)] <- "GL16_Input-WT"
treatment[grep("-Input-KO_", pnames)] <- "GL16_Input-KO"
treatment[grep("-Input-KO_", pnames)] <- "GL16_Input-KO"
treatment[grep("-MOZ-Asc-WT_", pnames)] <- "GL16_MOZ-Asc-WT"
treatment[grep("-MOZ-Asc-WT_", pnames)] <- "GL16_MOZ-Asc-WT"
treatment[grep("-MOZ-Asc-KO_", pnames)] <- "GL16_MOZ-Asc-KO"
treatment[grep("-MOZ-Asc-KO_", pnames)] <- "GL16_MOZ-Asc-KO"
treatment[grep("-H3K9-14-ac-WT_", pnames)] <- "GL16_H3K9-14-ac-WT"
treatment[grep("-H3K9-14-ac-WT_", pnames)] <- "GL16_H3K9-14-ac-WT"
treatment[grep("-H3K9-14-ac-KO_", pnames)] <- "GL16_H3K9-14-ac-KO"
treatment[grep("-H3K9-14-ac-KO_", pnames)] <- "GL16_H3K9-14-ac-KO"
treatment[grep("-H4k16ac-WT_", pnames)] <- "GL17_H4k16ac-WT"
treatment[grep("-H4k16ac-WT_", pnames)] <- "GL17_H4k16ac-WT"
treatment[grep("-H4k16ac-KO_", pnames)] <- "GL17_H4k16ac-KO"
treatment[grep("-H4k16ac-KO_", pnames)] <- "GL17_H4k16ac-KO"
treatment[grep("-H3k4me1-WT_", pnames)] <- "GL17_H3k4me1-WT"
treatment[grep("-H3k4me1-WT_", pnames)] <- "GL17_H3k4me1-WT"
treatment[grep("-H3k4me1-KO_", pnames)] <- "GL17_H3k4me1-KO"
treatment[grep("-H3k4me1-KO_", pnames)] <- "GL17_H3k4me1-KO"
treatment[grep("-Pol2-WT_", pnames)] <- "GL17_Pol2-WT"
treatment[grep("-Pol2-WT_", pnames)] <- "GL17_Pol2-WT"
treatment[grep("-Pol2-KO_", pnames)] <- "GL17_Pol2-KO"
treatment[grep("-Pol2-KO_", pnames)] <- "GL17_Pol2-KO"
treatment[grep("-H3k4me3-WT_", pnames)] <- "GL17_H3k4me3-WT"
treatment[grep("-H3k4me3-WT_", pnames)] <- "GL17_H3k4me3-WT"
treatment[grep("-H3k4me3-KO_", pnames)] <- "GL17_H3k4me3-KO"
treatment[grep("-H3k4me3-KO_", pnames)] <- "GL17_H3k4me3-KO"
treatment[grep("-H3k27me3-WT_", pnames)] <- "GL18_H3k27me3-WT"
treatment[grep("-H3k27me3-WT_", pnames)] <- "GL18_H3k27me3-WT"
treatment[grep("-H3k27me3-KO_", pnames)] <- "GL18_H3k27me3-KO"
treatment[grep("-H3k27me3-KO_", pnames)] <- "GL18_H3k27me3-KO"

groups <- factor(treatment, levels=unique(treatment))

sample.desc <- data.frame(filename=myfiles, key=groups)
write.table(sample.desc,file.path(output.dir,"sample.desc.txt"),quote=F,sep="\t",col.names=NA)
fac <- unique(as.character(sample.desc$key))


for(i in 1:length(fac)){
  ID = fac[i]
  cat("Merging ", ID, "\n", sep="")
  sel <- which(as.character(sample.desc[, "key"]) == ID)
  in.bam = paste(file.path(input.dir, as.character(sample.desc[sel, "filename"])), collapse=" ")  
  outpath.sampleName = file.path(dest.dir, ID)

#Why does the samtools command has the "-" symbol?? & samtools sort as well
  mergeBam.cmd <- paste("samtools merge -", in.bam, "| samtools sort -", outpath.sampleName, sep=" ")
  indexBam.cmd <- paste("samtools index ", outpath.sampleName, ".bam", sep="")

  system(mergeBam.cmd)
  system(indexBam.cmd)

  write(ID, file=command.file, sep="\n", append=TRUE)
  write(mergeBam.cmd, file=command.file, sep="\n", append=TRUE)
  write(indexBam.cmd, file=command.file, sep="\n", append=TRUE)
  write("\n", file=command.file, sep="\n", append=TRUE)
 
}
