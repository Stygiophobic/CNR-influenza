library(stringr)
#setwd("C:/Users/regueex/Desktop/test_merge_meta")
reName<-function(name){
  split_string<-strsplit(name,"_")
  sampleName<-split_string[[1]][1]
  return(sampleName)
}

reNameLoc<-function(name){
  split_string<-strsplit(name," / ")
  sampleName<-split_string[[1]][1]
  return(sampleName)
}

###########################################################################
argv <- commandArgs(TRUE)
raw_metadata<-read.csv2(as.character(argv[1]))

#raw_metadata<-read.csv2("metadata_raw.csv")

cleaned_metadata<-subset(raw_metadata,select=c("Seq_Id..HA.","Seq_Id..NA.",
                                               "Subtype","Lineage","Location","Host",
                                               "Collection_Date"))
colnames(cleaned_metadata)<-c("S4","S6","Subtype","Lineage","Location","Host","Date")

cleaned_metadata$Subtype<-paste0(cleaned_metadata$Subtype,cleaned_metadata$Lineage)

cleaned_metadata<-cleaned_metadata[,-4]
cleaned_metadata$Subtype<-str_replace(cleaned_metadata$Subtype, "H1N1pdm", "H1N1")
cleaned_metadata$Subtype<-str_replace(cleaned_metadata$Subtype, "BVictoria", "BVIC")
cleaned_metadata$Subtype<-str_replace(cleaned_metadata$Subtype, "BYamagata", "BYAM")

H1N1_metadata<-subset(cleaned_metadata,Subtype=="H1N1")
H3N2_metadata<-subset(cleaned_metadata,Subtype=="H3N2")
BVIC_metadata<-subset(cleaned_metadata,Subtype=="BVIC")
BYAM_metadata<-subset(cleaned_metadata,Subtype=="BYAM")

label<-c("strain","subtype","country","host","date")

H1N1_S4<-subset(H1N1_metadata,select=-c(S6))
H1N1_S6<-subset(H1N1_metadata,select=-c(S4))
colnames(H1N1_S4)<-label
colnames(H1N1_S6)<-label

H3N2_S4<-subset(H3N2_metadata,select=-c(S6))
H3N2_S6<-subset(H3N2_metadata,select=-c(S4))
colnames(H3N2_S4)<-label
colnames(H3N2_S6)<-label

BVIC_S4<-subset(BVIC_metadata,select=-c(S6))
BVIC_S6<-subset(BVIC_metadata,select=-c(S4))
colnames(BVIC_S4)<-label
colnames(BVIC_S6)<-label

BYAM_S4<-subset(BYAM_metadata,select=-c(S6))
BYAM_S6<-subset(BYAM_metadata,select=-c(S4))
colnames(BYAM_S4)<-label
colnames(BYAM_S6)<-label


###########################################################################
ref_H1N1<-read.csv2(as.character(argv[2]))
#ref_H1N1<-read.csv2("meta_H1N1.csv")

clean_data<-subset(ref_H1N1,select=c("HA.Segment_Id","NA.Segment_Id",
                                 "Location","Host","Collection_Date",
                                 "Isolate_Name",	"Subtype",	"Lineage"))
colnames(clean_data)<-c("H4","H6","Location","Host","Collection_Date",
                        "Isolate_Name",	"Subtype",	"Lineage")

clean_data$H4<-paste0(clean_data$Isolate_Name,"_HA")
clean_data$H6<-paste0(clean_data$Isolate_Name,"_NA")

clean_data<-clean_data[,-6]

clean_data$Subtype<-str_replace(clean_data$Subtype,"A / ","")
clean_data$Lineage<-str_replace(clean_data$Lineage,"pdm09","")
clean_data$Subtype<-paste0(clean_data$Subtype,clean_data$Lineage)
clean_data<-clean_data[,-7]

clean_data$Location<-sapply(as.character(clean_data$Location),reNameLoc)

clean_data$H4<-str_replace(clean_data$H4," ","_")
clean_data$H6<-str_replace(clean_data$H6," ","_")
label<-c("strain","subtype","country","host","date")

ref_H1N1_S4<-clean_data[,c(1,6,3,4,5)]
colnames(ref_H1N1_S4)<-label
ref_H1N1_S6<-clean_data[,c(2,6,3,4,5)]
colnames(ref_H1N1_S6)<-label


###########################################################################
ref_H3N2<-read.csv2(as.character(argv[3]))
#ref_H3N2<-read.csv2("meta_H3N2.csv")

clean_data<-subset(ref_H3N2,select=c("HA.Segment_Id","NA.Segment_Id",
                                     "Location","Host","Collection_Date",
                                     "Isolate_Name",	"Subtype",	"Lineage"))
colnames(clean_data)<-c("H4","H6","Location","Host","Collection_Date",
                        "Isolate_Name",	"Subtype",	"Lineage")

clean_data$H4<-paste0(clean_data$Isolate_Name,"_HA")
clean_data$H6<-paste0(clean_data$Isolate_Name,"_NA")

clean_data<-clean_data[,-6]

clean_data$Subtype<-str_replace(clean_data$Subtype,"A / ","")
clean_data<-clean_data[,-7]

clean_data$Location<-sapply(as.character(clean_data$Location),reNameLoc)

clean_data$H4<-str_replace(clean_data$H4," ","_")
clean_data$H6<-str_replace(clean_data$H6," ","_")
label<-c("strain","subtype","country","host","date")

ref_H3N2_S4<-clean_data[,c(1,6,3,4,5)]
colnames(ref_H3N2_S4)<-label
ref_H3N2_S6<-clean_data[,c(2,6,3,4,5)]
colnames(ref_H3N2_S6)<-label

###########################################################################
ref_BVIC<-read.csv2(as.character(argv[4]))
#ref_BVIC<-read.csv2("meta_BVIC.csv")

clean_data<-subset(ref_BVIC,select=c("HA.Segment_Id","NA.Segment_Id",
                                     "Location","Host","Collection_Date",
                                     "Isolate_Name",	"Subtype",	"Lineage"))
colnames(clean_data)<-c("H4","H6","Location","Host","Collection_Date",
                        "Isolate_Name",	"Subtype",	"Lineage")

clean_data$H4<-paste0(clean_data$Isolate_Name,"_HA")
clean_data$H6<-paste0(clean_data$Isolate_Name,"_NA")

clean_data<-clean_data[,-6]

clean_data$Subtype<-str_replace(clean_data$Subtype,"B","BVIC")
clean_data<-clean_data[,-7]

clean_data$Location<-sapply(as.character(clean_data$Location),reNameLoc)

clean_data$H4<-str_replace(clean_data$H4," ","_")
clean_data$H6<-str_replace(clean_data$H6," ","_")
label<-c("strain","subtype","country","host","date")

ref_BVIC_S4<-clean_data[,c(1,6,3,4,5)]
colnames(ref_BVIC_S4)<-label
ref_BVIC_S6<-clean_data[,c(2,6,3,4,5)]
colnames(ref_BVIC_S6)<-label

###########################################################################
ref_BYAM<-read.csv2(as.character(argv[5]))
#ref_BYAM<-read.csv2("meta_BYAM.csv")

clean_data<-subset(ref_BYAM,select=c("HA.Segment_Id","NA.Segment_Id",
                                     "Location","Host","Collection_Date",
                                     "Isolate_Name",	"Subtype",	"Lineage"))
colnames(clean_data)<-c("H4","H6","Location","Host","Collection_Date",
                        "Isolate_Name",	"Subtype",	"Lineage")

clean_data$H4<-paste0(clean_data$Isolate_Name,"_HA")
clean_data$H6<-paste0(clean_data$Isolate_Name,"_NA")

clean_data<-clean_data[,-6]

clean_data$Subtype<-str_replace(clean_data$Subtype,"B","BYAM")
clean_data<-clean_data[,-7]

clean_data$Location<-sapply(as.character(clean_data$Location),reNameLoc)

clean_data$H4<-str_replace(clean_data$H4," ","_")
clean_data$H6<-str_replace(clean_data$H6," ","_")
label<-c("strain","subtype","country","host","date")

ref_BYAM_S4<-clean_data[,c(1,6,3,4,5)]
colnames(ref_BYAM_S4)<-label
ref_BYAM_S6<-clean_data[,c(2,6,3,4,5)]
colnames(ref_BYAM_S6)<-label

############################################################################

H1N1_S4<-rbind(H1N1_S4,ref_H1N1_S4)
H1N1_S6<-rbind(H1N1_S6,ref_H1N1_S6)

H3N2_S4<-rbind(H3N2_S4,ref_H3N2_S4)
H3N2_S6<-rbind(H3N2_S6,ref_H3N2_S6)

BVIC_S4<-rbind(BVIC_S4,ref_BVIC_S4)
BVIC_S6<-rbind(BVIC_S6,ref_BVIC_S6)

BYAM_S4<-rbind(BYAM_S4,ref_BYAM_S4)
BYAM_S6<-rbind(BYAM_S6,ref_BYAM_S6)

############################################################################
write.table(H1N1_S4,file="temp_data/H1N1_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H1N1_S6,file="temp_data/H1N1_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H3N2_S4,file="temp_data/H3N2_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H3N2_S6,file="temp_data/H3N2_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BYAM_S4,file="temp_data/BYAM_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BYAM_S6,file="temp_data/BYAM_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BVIC_S4,file="temp_data/BVIC_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BVIC_S6,file="temp_data/BVIC_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
