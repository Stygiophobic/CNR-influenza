argv <- commandArgs(TRUE)
raw_metadata<-read.csv2(as.character(argv[1]))
#raw_metadata<-read.csv2("/home/hadrien/Bureau/CNR-influenza/temp_data/metadata_raw.csv")

cleaned_metadata<-subset(raw_metadata,select=c("Seq_Id..HA.","Seq_Id..NA.","Subtype","Lineage","Location","Host","Collection_Date"))
colnames(cleaned_metadata)<-c("S4","S6","subtype","lineage","country","host","date")

tableS4<-cleaned_metadata[,c(1,3,4,5,6,7)]
colnames(tableS4)[1]<-"strain"
tableS4$segment<-"S4"

tableS6<-cleaned_metadata[,c(2,3,4,5,6,7)]
colnames(tableS6)[1]<-"strain"
tableS6$segment<-"S6"

full_data<-rbind(tableS4,tableS6)
full_data$subtype<-paste0(full_data$subtype,full_data$lineage)
full_data<-full_data[,-3]

write.table(full_data,file="temp_data/metadata.tsv",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
