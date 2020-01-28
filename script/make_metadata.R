argv <- commandArgs(TRUE)

#raw_metadata<-snakemake@input[[1]]
#H1N1_S4_output<-snakemake@output[[1]]
#H1N1_S6_output<-snakemake@output[[2]]
#H3N2_S4_output<-snakemake@output[[3]]
#H3N2_S4_output<-snakemake@output[[4]]
#B_S4_output<-snakemake@output[[5]]
#B_S4_output<-snakemake@output[[6]]
#snakemake@config[["myparam"]]

raw_metadata<-read.csv2(as.character(argv[1]))
#H1N1_S4_output<-as.character(argv[2])
#H1N1_S6_output<-as.character(argv[3])
#H3N2_S4_output<-as.character(argv[4])
#H3N2_S6_output<-as.character(argv[5])
#B_S4_output<-as.character(argv[6])
#B_S6_output<-as.character(argv[7])

#raw_metadata<-read.csv2("/home/hadrien/Bureau/CNR-influenza/temp_data/metadata_raw.csv")

cleaned_metadata<-subset(raw_metadata,select=c("Seq_Id..HA.","Seq_Id..NA.","Subtype","Location","Host","Collection_Date"))
colnames(cleaned_metadata)<-c("S4","S6","Subtype","Location","Host","Date")

H1N1_metadata<-subset(cleaned_metadata,Subtype=="H1N1")
H3N2_metadata<-subset(cleaned_metadata,Subtype=="H3N2")
B_metadata<-subset(cleaned_metadata,Subtype=="B")

label<-c("strain","country","host","date")

H1N1_S4<-subset(H1N1_metadata,select=-c(S6,Subtype))
H1N1_S6<-subset(H1N1_metadata,select=-c(S4,Subtype))
colnames(H1N1_S4)<-label
colnames(H1N1_S6)<-label

H3N2_S4<-subset(H3N2_metadata,select=-c(S6,Subtype))
H3N2_S6<-subset(H3N2_metadata,select=-c(S4,Subtype))
colnames(H3N2_S4)<-label
colnames(H3N2_S6)<-label

B_S4<-subset(B_metadata,select=-c(S6,Subtype))
B_S6<-subset(B_metadata,select=-c(S4,Subtype))
colnames(B_S4)<-label
colnames(B_S6)<-label

#write.table(H1N1_S4,file=H1N1_S4,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#write.table(H1N1_S6,file=H1N1_S6,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#write.table(H3N2_S4,file=H3N2_S4,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#write.table(H3N2_S6,file=H3N2_S6,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#write.table(B_S4,file=B_S4,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
#write.table(B_S6,file=B_S6,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")

write.table(H1N1_S4,file="temp_data/H1N1S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H1N1_S6,file="temp_data/H1N1S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H3N2_S4,file="temp_data/H3N2S4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(H3N2_S6,file="temp_data/H3N2S6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(B_S4,file="temp_data/BS4.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
write.table(B_S6,file="temp_data/BS6.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
