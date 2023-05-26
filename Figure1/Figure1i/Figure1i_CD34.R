#Figure1i
#Author: Alba Rodriguez-Meira
##############################

#Summarize karyotipic abnormalities in TP53-mutant cells from IF0131 and GR001 at baseline (in CD34+ cells)
#This data is obtained from the inferCNV readouts.

IF0131<-read.table(file="CD34_inferCNV/cnv_summary_IF0131_subclusters_genotype.txt",sep = "\t",header =T)
table(IF0131$cluster_corrected)
#chr3_chr5
584/length(IF0131$cluster_corrected) #0.5402405
#chr3_chr5_chr7
248/length(IF0131$cluster_corrected) #0.2294172
####################################
GR001<-read.table(file="CD34_inferCNV/cnv_summary_GR001_subclusters_genotype.txt",sep = "\t",header =T)
table(GR001$cluster_corrected)
#chr5_chr17_chr21
753/length(GR001$cluster_corrected) #0.7440711
#chr5_chr7_chr9_chr21
113/length(GR001$cluster_corrected) #0.1116601
