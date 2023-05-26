#FigureS3a - quantification of index sorting data; Revision 1

################################################
metadata<-read.table('../../../Sources/metadata_MPNAMLp53_with_index_genotype.revised.txt',header = T,sep = '\t')
################################################

table(metadata$stage)

#MEPs (CD38>value,CD45RA<value,CD123<value)

GR003.MEP<-subset(metadata,stage=="GR003_AP" & CD38>1500 & CD45RA<3000 & CD123<800)
IF0131.MEP<-subset(metadata,stage=="IF0131" & CD38>10000 & CD45RA<7000 & CD123<1000)
SB5702.MEP<-subset(metadata,stage=="SB5702" & CD38>10000 & CD45RA<7000 & CD123<1000)
GR001.MEP<-subset(metadata,stage=="GR001" & CD38>5000 & CD45RA<5000 & CD123<1000)
GR005.MEP<-subset(metadata,stage=="GR005_AP" & CD38>370 & CD45RA<3500 & CD123<1000)
GR007.MEP<-subset(metadata,stage=="GR007_AP" & CD38>2000 & CD45RA<4000 & CD123<1000)
GR006.MEP<-subset(metadata,stage=="GR006_AP" & CD38>1500 & CD45RA<5000 & CD123<1200)
PM4908.MEP<-subset(metadata,stage=="PM4908" & CD38>10000 & CD45RA<5000 & CD123<1000)
GR002.MEP<-subset(metadata,stage=="GR002" & CD38>1400 & CD45RA<3000 & CD123<230)
IF0318.MEP<-subset(metadata,stage=="IF0318" & CD38>6000 & CD45RA<5000 & CD123<1000)
DT5210.MEP<-subset(metadata,stage=="DT5210" & CD38>5000 & CD45RA<4000 & CD123<800)
GST010.MEP<-subset(metadata,stage=="GST010" & CD38>5000 & CD45RA<4000 & CD123<800)
JB4211.MEP<-subset(metadata,stage=="JB4211" & CD38>5000 & CD45RA<4000 & CD123<800)
IF0393.MEP<-subset(metadata,stage=="IF0393" & CD38>6000 & CD45RA<5000 & CD123<1000)
GR004_AP.MEP<-subset(metadata,stage=="GR004_AP" & CD38>1500 & CD45RA<3000 & CD123<800)
GH001_003.MEP<-subset(metadata,stage=="GH001_003" & CD38>1500 & CD45RA<3000 & CD123<800)

MEPs<-c(as.character(GR003.MEP$genotype.classification),
        as.character(IF0131.MEP$genotype.classification),
        as.character(SB5702.MEP$genotype.classification),
        as.character(GR001.MEP$genotype.classification),
        as.character(GR005.MEP$genotype.classification),
        as.character(GR007.MEP$genotype.classification),
        as.character(GR006.MEP$genotype.classification),
        as.character(PM4908.MEP$genotype.classification),
        as.character(GR002.MEP$genotype.classification),
        as.character(IF0318.MEP$genotype.classification),
        as.character(DT5210.MEP$genotype.classification),
        as.character(GST010.MEP$genotype.classification),
        as.character(JB4211.MEP$genotype.classification),
        as.character(IF0393.MEP$genotype.classification),
        as.character(GR004_AP.MEP$genotype.classification),
        as.character(GH001_003.MEP$genotype.classification))

#GMPs (CD38>value,CD45RA>value,CD123>value)

GR003.GMP<-subset(metadata,stage=="GR003_AP" & CD38>1500 & CD45RA>3000 & CD123>800)
IF0131.GMP<-subset(metadata,stage=="IF0131" & CD38>10000 & CD45RA>7000 & CD123>1000)
SB5702.GMP<-subset(metadata,stage=="SB5702" & CD38>10000 & CD45RA>7000 & CD123>1000)
GR001.GMP<-subset(metadata,stage=="GR001" & CD38>5000 & CD45RA>5000 & CD123>1000)
GR005.GMP<-subset(metadata,stage=="GR005_AP" & CD38>370 & CD45RA>3500 & CD123>1000)
GR007.GMP<-subset(metadata,stage=="GR007_AP" & CD38>2000 & CD45RA>4000 & CD123>1000)
GR006.GMP<-subset(metadata,stage=="GR006_AP" & CD38>1500 & CD45RA>5000 & CD123>1200)
PM4908.GMP<-subset(metadata,stage=="PM4908" & CD38>10000 & CD45RA>5000 & CD123>1000)
GR002.GMP<-subset(metadata,stage=="GR002" & CD38>1400 & CD45RA>3000 & CD123>230)
IF0318.GMP<-subset(metadata,stage=="IF0318" & CD38>6000 & CD45RA>5000 & CD123>1000)
DT5210.GMP<-subset(metadata,stage=="DT5210" & CD38>5000 & CD45RA>4000 & CD123>800)
GST010.GMP<-subset(metadata,stage=="GST010" & CD38>5000 & CD45RA>4000 & CD123>800)
JB4211.GMP<-subset(metadata,stage=="JB4211" & CD38>5000 & CD45RA>4000 & CD123>800)
IF0393.GMP<-subset(metadata,stage=="IF0393" & CD38>6000 & CD45RA>5000 & CD123>1000)
GR004_AP.GMP<-subset(metadata,stage=="GR004_AP" & CD38>1500 & CD45RA>3000 & CD123>800)
GH001_003.GMP<-subset(metadata,stage=="GH001_003" & CD38>1500 & CD45RA>3000 & CD123>800)

GMPs<-c(as.character(GR003.GMP$genotype.classification),
        as.character(IF0131.GMP$genotype.classification),
        as.character(SB5702.GMP$genotype.classification),
        as.character(GR001.GMP$genotype.classification),
        as.character(GR005.GMP$genotype.classification),
        as.character(GR007.GMP$genotype.classification),
        as.character(GR006.GMP$genotype.classification),
        as.character(PM4908.GMP$genotype.classification),
        as.character(GR002.GMP$genotype.classification),
        as.character(IF0318.GMP$genotype.classification),
        as.character(DT5210.GMP$genotype.classification),
        as.character(GST010.GMP$genotype.classification),
        as.character(JB4211.GMP$genotype.classification),
        as.character(IF0393.GMP$genotype.classification),
        as.character(GR004_AP.GMP$genotype.classification),
        as.character(GH001_003.GMP$genotype.classification))

#CMPs (CD38>value,CD45RA<value,CD123>value)

GR003.CMP<-subset(metadata,stage=="GR003_AP" & CD38>1500 & CD45RA<3000 & CD123>800)
IF0131.CMP<-subset(metadata,stage=="IF0131" & CD38>10000 & CD45RA<7000 & CD123>1000)
SB5702.CMP<-subset(metadata,stage=="SB5702" & CD38>10000 & CD45RA<7000 & CD123>1000)
GR001.CMP<-subset(metadata,stage=="GR001" & CD38>5000 & CD45RA<5000 & CD123>1000)
GR005.CMP<-subset(metadata,stage=="GR005_AP" & CD38>370 & CD45RA<3500 & CD123>1000)
GR007.CMP<-subset(metadata,stage=="GR007_AP" & CD38>2000 & CD45RA<4000 & CD123>1000)
GR006.CMP<-subset(metadata,stage=="GR006_AP" & CD38>1500 & CD45RA<5000 & CD123>1200)
PM4908.CMP<-subset(metadata,stage=="PM4908" & CD38>10000 & CD45RA<5000 & CD123>1000)
GR002.CMP<-subset(metadata,stage=="GR002" & CD38>1400 & CD45RA<3000 & CD123>230)
IF0318.CMP<-subset(metadata,stage=="IF0318" & CD38>6000 & CD45RA<5000 & CD123>1000)
DT5210.CMP<-subset(metadata,stage=="DT5210" & CD38>5000 & CD45RA<4000 & CD123>800)
GST010.CMP<-subset(metadata,stage=="GST010" & CD38>5000 & CD45RA<4000 & CD123>800)
JB4211.CMP<-subset(metadata,stage=="JB4211" & CD38>5000 & CD45RA<4000 & CD123>800)
IF0393.CMP<-subset(metadata,stage=="IF0393" & CD38>6000 & CD45RA<5000 & CD123>1000)
GR004_AP.CMP<-subset(metadata,stage=="GR004_AP" & CD38>1500 & CD45RA<3000 & CD123>800)
GH001_003.CMP<-subset(metadata,stage=="GH001_003" & CD38>1500 & CD45RA<3000 & CD123>800)

CMPs<-c(as.character(GR003.CMP$genotype.classification),
        as.character(IF0131.CMP$genotype.classification),
        as.character(SB5702.CMP$genotype.classification),
        as.character(GR001.CMP$genotype.classification),
        as.character(GR005.CMP$genotype.classification),
        as.character(GR007.CMP$genotype.classification),
        as.character(GR006.CMP$genotype.classification),
        as.character(PM4908.CMP$genotype.classification),
        as.character(GR002.CMP$genotype.classification),
        as.character(IF0318.CMP$genotype.classification),
        as.character(DT5210.CMP$genotype.classification),
        as.character(GST010.CMP$genotype.classification),
        as.character(JB4211.CMP$genotype.classification),
        as.character(IF0393.CMP$genotype.classification),
        as.character(GR004_AP.CMP$genotype.classification),
        as.character(GH001_003.CMP$genotype.classification))

#LMPPs (CD38<value,CD45RA>value,CD90<value)

GR003.LMPP<-subset(metadata,stage=="GR003_AP" & CD38<1500 & CD45RA>3000 & CD90<1500)
IF0131.LMPP<-subset(metadata,stage=="IF0131" & CD38<10000 & CD45RA>7000 & CD90<10000)
SB5702.LMPP<-subset(metadata,stage=="SB5702" & CD38<10000 & CD45RA>7000 & CD90<10000)
GR001.LMPP<-subset(metadata,stage=="GR001" & CD38<5000 & CD45RA>5000 & CD90<10000)
GR005.LMPP<-subset(metadata,stage=="GR005_AP" & CD38<370 & CD45RA>3500 & CD90<10000)
GR007.LMPP<-subset(metadata,stage=="GR007_AP" & CD38<2000 & CD45RA>4000 & CD90<1500)
GR006.LMPP<-subset(metadata,stage=="GR006_AP" & CD38<1500 & CD45RA>5000 & CD90<10000)
PM4908.LMPP<-subset(metadata,stage=="PM4908" & CD38<10000 & CD45RA>5000 & CD90<10000)
GR002.LMPP<-subset(metadata,stage=="GR002" & CD38<1400 & CD45RA>3000 & CD90<10000)
IF0318.LMPP<-subset(metadata,stage=="IF0318" & CD38<6000 & CD45RA>5000 & CD90<5000)
DT5210.LMPP<-subset(metadata,stage=="DT5210" & CD38<5000 & CD45RA>4000 & CD90<10000)
GST010.LMPP<-subset(metadata,stage=="GST010" & CD38<5000 & CD45RA>4000 & CD90<10000)
JB4211.LMPP<-subset(metadata,stage=="JB4211" & CD38<5000 & CD45RA>4000 & CD90<10000)
IF0393.LMPP<-subset(metadata,stage=="IF0393" & CD38<6000 & CD45RA>5000 & CD90<5000)
GR004_AP.LMPP<-subset(metadata,stage=="GR004_AP" & CD38<1500 & CD45RA>3000 & CD90<1500)
GH001_003.LMPP<-subset(metadata,stage=="GH001_003" & CD38<1500 & CD45RA>3000 & CD90<1500)

LMPPs<-c(as.character(GR003.LMPP$genotype.classification),
         as.character(IF0131.LMPP$genotype.classification),
         as.character(SB5702.LMPP$genotype.classification),
         as.character(GR001.LMPP$genotype.classification),
         as.character(GR005.LMPP$genotype.classification),
         as.character(GR007.LMPP$genotype.classification),
         as.character(GR006.LMPP$genotype.classification),
         as.character(PM4908.LMPP$genotype.classification),
         as.character(GR002.LMPP$genotype.classification),
         as.character(IF0318.LMPP$genotype.classification),
         as.character(DT5210.LMPP$genotype.classification),
         as.character(GST010.LMPP$genotype.classification),
         as.character(JB4211.LMPP$genotype.classification),
         as.character(IF0393.LMPP$genotype.classification),
         as.character(GR004_AP.LMPP$genotype.classification),
         as.character(GH001_003.LMPP$genotype.classification))

#MPPs (CD38<value,CD45RA<value,CD90<value)

GR003.MPP<-subset(metadata,stage=="GR003_AP" & CD38<1500 & CD45RA<3000 & CD90<1500)
IF0131.MPP<-subset(metadata,stage=="IF0131" & CD38<10000 & CD45RA<7000 & CD90<10000)
SB5702.MPP<-subset(metadata,stage=="SB5702" & CD38<10000 & CD45RA<7000 & CD90<10000)
GR001.MPP<-subset(metadata,stage=="GR001" & CD38<5000 & CD45RA<5000 & CD90<10000)
GR005.MPP<-subset(metadata,stage=="GR005_AP" & CD38<370 & CD45RA<3500 & CD90<10000)
GR007.MPP<-subset(metadata,stage=="GR007_AP" & CD38<2000 & CD45RA<4000 & CD90<1500)
GR006.MPP<-subset(metadata,stage=="GR006_AP" & CD38<1500 & CD45RA<5000 & CD90<10000)
PM4908.MPP<-subset(metadata,stage=="PM4908" & CD38<10000 & CD45RA<5000 & CD90<10000)
GR002.MPP<-subset(metadata,stage=="GR002" & CD38<1400 & CD45RA<3000 & CD90<10000)
IF0318.MPP<-subset(metadata,stage=="IF0318" & CD38<6000 & CD45RA<5000 & CD90<5000)
DT5210.MPP<-subset(metadata,stage=="DT5210" & CD38<5000 & CD45RA<4000 & CD90<10000)
GST010.MPP<-subset(metadata,stage=="GST010" & CD38<5000 & CD45RA<4000 & CD90<10000)
JB4211.MPP<-subset(metadata,stage=="JB4211" & CD38<5000 & CD45RA<4000 & CD90<10000)
IF0393.MPP<-subset(metadata,stage=="IF0393" & CD38<6000 & CD45RA<5000 & CD90<5000)
GR004_AP.MPP<-subset(metadata,stage=="GR004_AP" & CD38<1500 & CD45RA<3000 & CD90<1500)
GH001_003.MPP<-subset(metadata,stage=="GH001_003" & CD38<1500 & CD45RA<3000 & CD90<1500)

MPPs<-c(as.character(GR003.MPP$genotype.classification),
        as.character(IF0131.MPP$genotype.classification),
        as.character(SB5702.MPP$genotype.classification),
        as.character(GR001.MPP$genotype.classification),
        as.character(GR005.MPP$genotype.classification),
        as.character(GR007.MPP$genotype.classification),
        as.character(GR006.MPP$genotype.classification),
        as.character(PM4908.MPP$genotype.classification),
        as.character(GR002.MPP$genotype.classification),
        as.character(IF0318.MPP$genotype.classification),
        as.character(DT5210.MPP$genotype.classification),
        as.character(GST010.MPP$genotype.classification),
        as.character(JB4211.MPP$genotype.classification),
        as.character(IF0393.MPP$genotype.classification),
        as.character(GR004_AP.MPP$genotype.classification),
        as.character(GH001_003.MPP$genotype.classification))

#HSCs (CD38<value,CD45RA<value,CD90>value)

GR003.HSC<-subset(metadata,stage=="GR003_AP" & CD38<1500 & CD45RA<3000 & CD90>1500)
IF0131.HSC<-subset(metadata,stage=="IF0131" & CD38<10000 & CD45RA<7000 & CD90>10000)
SB5702.HSC<-subset(metadata,stage=="SB5702" & CD38<10000 & CD45RA<7000 & CD90>10000)
GR001.HSC<-subset(metadata,stage=="GR001" & CD38<5000 & CD45RA<5000 & CD90>10000)
GR005.HSC<-subset(metadata,stage=="GR005_AP" & CD38<370 & CD45RA<3500 & CD90>10000)
GR007.HSC<-subset(metadata,stage=="GR007_AP" & CD38<2000 & CD45RA<4000 & CD90>1500)
GR006.HSC<-subset(metadata,stage=="GR006_AP" & CD38<1500 & CD45RA<5000 & CD90>10000)
PM4908.HSC<-subset(metadata,stage=="PM4908" & CD38<10000 & CD45RA<5000 & CD90>10000)
GR002.HSC<-subset(metadata,stage=="GR002" & CD38<1400 & CD45RA<3000 & CD90>10000)
IF0318.HSC<-subset(metadata,stage=="IF0318" & CD38<6000 & CD45RA<5000 & CD90>5000)
DT5210.HSC<-subset(metadata,stage=="DT5210" & CD38<5000 & CD45RA<4000 & CD90>10000)
GST010.HSC<-subset(metadata,stage=="GST010" & CD38<5000 & CD45RA<4000 & CD90>10000)
JB4211.HSC<-subset(metadata,stage=="JB4211" & CD38<5000 & CD45RA<4000 & CD90>10000)
IF0393.HSC<-subset(metadata,stage=="IF0393" & CD38<6000 & CD45RA<5000 & CD90>5000)
GR004_AP.HSC<-subset(metadata,stage=="GR004_AP" & CD38<1500 & CD45RA<3000 & CD90>1500)
GH001_003.HSC<-subset(metadata,stage=="GH001_003" & CD38<1500 & CD45RA<3000 & CD90>1500)

HSCs<-c(as.character(GR003.HSC$genotype.classification),
        as.character(IF0131.HSC$genotype.classification),
        as.character(SB5702.HSC$genotype.classification),
        as.character(GR001.HSC$genotype.classification),
        as.character(GR005.HSC$genotype.classification),
        as.character(GR007.HSC$genotype.classification),
        as.character(GR006.HSC$genotype.classification),
        as.character(PM4908.HSC$genotype.classification),
        as.character(GR002.HSC$genotype.classification),
        as.character(IF0318.HSC$genotype.classification),
        as.character(DT5210.HSC$genotype.classification),
        as.character(GST010.HSC$genotype.classification),
        as.character(JB4211.HSC$genotype.classification),
        as.character(IF0393.HSC$genotype.classification),
        as.character(GR004_AP.HSC$genotype.classification),
        as.character(GH001_003.HSC$genotype.classification))

table(HSCs)
table(MPPs)
table(LMPPs)
table(CMPs)
table(MEPs)
table(GMPs)

