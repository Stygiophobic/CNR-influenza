argv <- commandArgs(TRUE)
raw_metadata<-read.csv2(as.character(argv[1]))
#subtype<-as.character(argv[2])
#setwd("C:/Users/regueex/Desktop")
#raw_metadata<-read.csv2("C:/Users/regueex/Desktop/metadata_raw.csv")

cleaned_metadata<-subset(raw_metadata,select=c("Seq_Id..HA.","Seq_Id..NA.",
                                               "Subtype","Lineage","Location","Host",
                                               "Collection_Date"))
colnames(cleaned_metadata)<-c("S4","S6","Subtype","Lineage","Location","Host","Date")

cleaned_metadata$Subtype<-paste0(cleaned_metadata$Subtype,cleaned_metadata$Lineage)

cleaned_metadata<-cleaned_metadata[,-4]

H1N1_metadata<-subset(cleaned_metadata,Subtype=="H1N1pdm")
H3N2_metadata<-subset(cleaned_metadata,Subtype=="H3N2")
BVIC_metadata<-subset(cleaned_metadata,Subtype=="BVictoria")
BYAM_metadata<-subset(cleaned_metadata,Subtype=="BYamagata")

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

write.table(H1N1_S4,file="temp_data/H1N1_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H1N1_S6,file="temp_data/H1N1_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H3N2_S4,file="temp_data/H3N2_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H3N2_S6,file="temp_data/H3N2_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BYAM_S4,file="temp_data/BYAM_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BYAM_S6,file="temp_data/BYAM_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BVIC_S4,file="temp_data/BVIC_S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(BVIC_S6,file="temp_data/BVIC_S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")


#if (subtype=="H1N1S4") write.table(H1N1_S4,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#if (subtype=="H1N1S6") write.table(H1N1_S6,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#if (subtype=="H3N2S4") write.table(H3N2_S4,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#if (subtype=="H3N2S6") write.table(H3N2_S6,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#if (subtype=="BYAMS4") write.table(BYAM_S4,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#if (subtype=="BYAMS6") write.table(BYAM_S6,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#if (subtype=="BVICS4") write.table(BVIC_S4,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#if (subtype=="BVICS6") write.table(BVIC_S6,file=paste0(subtype,".tsv"),col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
