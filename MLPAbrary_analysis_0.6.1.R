#!/usr/bin/Rscript
###MLPAbrary_analysis.R###
#Takes in genotyping information produced by AmpliVar and creates a summary report for CNV

VERSION="MLPAbrary_analysis.R Version: 0.6.1 Implemented by Olga Kondrashova 26/05/2015. Updated on 20/11/2015"

#DEFINE Functions
moveMe <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}


######MAIN########
suppressPackageStartupMessages(library('ggplot2'));
suppressPackageStartupMessages(library('getopt'));
suppressPackageStartupMessages(library('matrixStats'))

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).

optmx <- matrix(c(
  'help', 'h', 0, "logical", "Print this helpful helpscreen",
  'InputDir', 'I', 1, "character", "Required: input directory of AmpliVar analysis with genotype files",
  'ControlDir','C', 1, "character", "Optional: directory with control samples - genotype files",
  'OutputDir', 'O', 1, "character", "Required: output directory for summary reports",
  'Suspects', 'S', 1, "character", "Required: usual suspects file used for MLPAbrary analysis",
  'verbose', 'v', 0, "logical", "Print progress messages",
  'Version', 'V',0, "logical","Display script Version"
),ncol=5,byrow=T)

opt = getopt(optmx);

#usage statement
usage <- function(){
  getopt(optmx,command="MLPAbrary_analysis.R",usage=T);
}


# #Some example options for debugging
# opt=list()
# opt$InputDir="~/data/miseq/allocate/mlpabry/testing/"
# opt$OutputDir="~/data/miseq/allocate/mlpabry/testing_output/"
# opt$ControlDir="~/data/miseq/allocate/mlpabry/testing_controls/"
# opt$Suspects="~/data/resources/mlpabry/MLPAbry_middle_long_suspects.txt"
# opt$verbose=T
# opt$print=T

## Some error handling and option setting
if ( !is.null(opt$help) ){ cat("\n",usage(),"\n"); q(status=0)};
if ( !is.null(opt$Version) ){ cat(VERSION,"\n"); q(status=0)};
if ( is.null(opt$InputDir) ){ cat("\nError: No Input directory defined \n"); cat(usage(), "\n"); q(status=1)};
if ( is.null(opt$OutputDir) ){ cat("\nError: No Output directory defined \n"); cat(usage(), "\n"); q(status=1)};
if ( is.null(opt$Suspects) ){ cat("\nError: No suspects file defined \n"); cat(usage(), "\n"); q(status=1)};
if (!is.null(opt$verbose)) {cat(paste(Sys.time(),"Starting MLPAbrary analysis for",opt$InputDir,"\n"))
                            cat(VERSION,"\n")
                            cat("Options:\n")
                            print(opt)
}

####load the MLPAbrary data

#list files in the directory
if ( is.null(opt$ControlDir) ) {
  filenames=list.files(path=opt$InputDir,full.names=FALSE,recursive=FALSE)
} else {
  filenames=c(list.files(path=opt$InputDir,full.names=FALSE,recursive=FALSE),
              list.files(path=opt$Control,full.names=FALSE,recursive=FALSE))
}


#read in label columns
mlpabrary=read.delim(opt$Suspects,header=FALSE)[c(1,3)]
colnames(mlpabrary)=c("probe","chr_pos")

#create dataframe for mlpabrary results
for (i in 1:length(filenames)) { 
  name<-filenames[i]
  name2=sub("_merged_seqprep.fastq.fna_1_num_grp_Genotypes.txt","_cov",name)
  if (is.element(name,list.files(path=opt$InputDir,full.names=FALSE,recursive=FALSE)) ) {
    mlpabrary[,paste(name2)]<-read.delim(paste(opt$InputDir,"/",name,sep=""))[2]
  } else {
    mlpabrary[,paste(name2)]<-read.delim(paste(opt$ControlDir,"/",name,sep=""))[2]
  }
  mlpabrary[,paste(name2)]=as.numeric(mlpabrary[,paste(name2)])
}

#create subsets for mix A and mix B
mlpabrary_MixA=data.frame(mlpabrary[grepl("mixa",mlpabrary$probe,ignore.case=TRUE),
                                    grepl("mixa|probe|chr_pos",colnames(mlpabrary),ignore.case=TRUE)],
                          row.names=NULL,check.names=FALSE)
mlpabrary_MixB=data.frame(mlpabrary[grepl("mixb",mlpabrary$probe,ignore.case=TRUE),
                                    grepl("mixb|probe|chr_pos",colnames(mlpabrary),ignore.case=TRUE)],
                          row.names=NULL,check.names=FALSE)

#create normalised mlpabrary results for mix A and mix B
for (i in 3:length(colnames(mlpabrary_MixA))) { 
  name=colnames(mlpabrary_MixA)[i]
  name2=sub("_cov","_norm",name)
  control_probes=mean(mlpabrary_MixA[grep(".*ctrl$",mlpabrary_MixA$probe,ignore.case=T),name])
  mlpabrary_MixA[,name2]=mlpabrary_MixA[,name]/control_probes
}
for (i in 3:length(colnames(mlpabrary_MixB))) { 
  name=colnames(mlpabrary_MixB)[i]
  name2=sub("_cov","_norm",name)
  control_probes=mean(mlpabrary_MixB[grep(".*ctrl$",mlpabrary_MixB$probe,ignore.case=T),name])
  mlpabrary_MixB[,name2]=mlpabrary_MixB[,name]/control_probes
}

#Calculate ratios for mixA and mixB samples
for (i in 1:length(colnames(mlpabrary_MixA))) { 
  if (grepl("_norm$",colnames(mlpabrary_MixA)[i],ignore.case=T)==TRUE){
    name=colnames(mlpabrary_MixA)[i]
    if (length(grep("^.*control.*norm$",colnames(mlpabrary_MixA),ignore.case=TRUE))>1){
      mlpabrary_MixA$control_mean_mixA=rowMeans(mlpabrary_MixA[,grep("^.*control.*norm$",colnames(mlpabrary_MixA),ignore.case=T)])
      mlpabrary_MixA[,paste(name,"_ratio",sep="")]=(mlpabrary_MixA[,i]/mlpabrary_MixA$control_mean_mixA)
    }else{
      mlpabrary_MixA$control_mean_mixA=(mlpabrary_MixA[,grep("^.*control.*norm$",colnames(mlpabrary_MixA),
                                                             ignore.case=TRUE)]) 
      mlpabrary_MixA[,paste(name,"_ratio",sep="")]=
        mlpabrary_MixA[,i]/mlpabrary_MixA$control_mean_mixA 
    }
  }
}
for (i in 1:length(colnames(mlpabrary_MixB))) { 
  if (grepl("_norm$",colnames(mlpabrary_MixB)[i],ignore.case=T)==TRUE){
    name=colnames(mlpabrary_MixB)[i]
    if (length(grep("^.*control.*norm$",colnames(mlpabrary_MixB),ignore.case=TRUE))>1){
      mlpabrary_MixB$control_mean_MixB=rowMeans(mlpabrary_MixB[,grep("^.*control.*norm$",colnames(mlpabrary_MixB),ignore.case=T)])
      mlpabrary_MixB[,paste(name,"_ratio",sep="")]=(mlpabrary_MixB[,i]/mlpabrary_MixB$control_mean_MixB)
    }else{
      mlpabrary_MixB$control_mean_MixB=(mlpabrary_MixB[,grep("^.*control.*norm$",colnames(mlpabrary_MixB),
                                                             ignore.case=TRUE)])
      mlpabrary_MixB[,paste(name,"_ratio",sep="")]=
        mlpabrary_MixB[,i]/mlpabrary_MixB$control_mean_MixB
    }
  }
}


#merge back to one dataframe
final_values=merge(mlpabrary_MixA,mlpabrary_MixB, by=c("probe","chr_pos"),all = TRUE)

#Mean NTC cov
for (i in 1:length(colnames(final_values))){
  if (isTRUE(grepl("^.*NTC.*MixA.*cov",colnames(final_values)[i],ignore.case=TRUE))) {
    mean_cov_NTC_mixA=mean(final_values[,i],na.rm=TRUE)
  } else if (isTRUE(grepl("^.*NTC.*MixB.*cov",colnames(final_values)[i],ignore.case=TRUE)))  {
    mean_cov_NTC_mixB=mean(final_values[,i],na.rm=TRUE)    
  }
}  

#Get mean for two probe sets and difference
for (i in 1:length(colnames(final_values))){
  if (isTRUE(grepl("^.*mixA.*_norm_ratio$",colnames(final_values)[i],ignore.case=TRUE))) {
    samplename_report_mixA=sub("_norm_ratio","",colnames(final_values)[i])
    meancov_mixA_name=sub("_norm_ratio","_cov",colnames(final_values)[i])
    samplename_mixA=unlist(strsplit(colnames(final_values)[i],"[_]"))[5]
#   Little fix up for case specific differences
    samplename_mixA=sub("mixA","MixA",samplename_mixA)
    samplename_mixB=sub("MixA","MixB",samplename_mixA)
    samplename=unlist(strsplit(samplename_mixA,"[-]"))[1]
    grepl_expr_mixB=paste("^.*",samplename_mixB,".*ratio",sep="")
    for (n in 1:length(colnames(final_values))) {
      if(isTRUE(grepl(grepl_expr_mixB,colnames(final_values)[n],ignore.case=TRUE))){
        meancov_mixB_name=sub("_norm_ratio","_cov",colnames(final_values)[n])
        samplename_report_mixB=sub("_norm_ratio","",colnames(final_values)[n])
        final_values[,paste(samplename,"_meanratio",sep="")]=rowMeans(final_values[,c(i,n)], na.rm = TRUE)
        final_values[,paste(samplename,"_stdev",sep="")]=rowSds(data.matrix(final_values[,c(i,n)],rownames.force = NA),na.rm=TRUE)

      }
    }
  }
}

# Combine different probes for each exon (BRCA1/2 and PTEN)
final_values$exon=apply(sapply(strsplit(as.vector(final_values$probe), "[_]"), `[`, c(2,3)),2,paste,collapse="_") 
combined_final_values=data.frame(exon=unique(final_values$exon))

for (i in 1:length(colnames(final_values))){
  if (isTRUE(grepl("*_meanratio",colnames(final_values)[i],ignore.case=TRUE))) {
    name=sub("_meanratio","_combinedratio",colnames(final_values)[i])
    name2=sub("_meanratio","_stdev_exon",colnames(final_values)[i])
    c=aggregate(final_values[,i],by=list(exon=final_values$exon),function(x) c(mean=mean(x),stdev=sd(x)))
    combined_final_values[,c("exon",name,name2)] = c(c$exon, c$x[,"mean"], c$x[,"stdev"])
    #round up the output values  
    combined_final_values[,name]=round(as.numeric(combined_final_values[,name]),digits=2) 
    combined_final_values[,name2]=round(as.numeric(combined_final_values[,name2]),digits=2) 
  }
  if (isTRUE(grepl("*_stdev",colnames(final_values)[i],ignore.case=TRUE))) {
    name=sub("_stdev","_stdev_mix",colnames(final_values)[i])
    combined_final_values[,c("exon",name)]=aggregate(final_values[,i],by=list(exon=final_values$exon),FUN=mean)
    #round up the output values  
    combined_final_values[,name]=round(as.numeric(combined_final_values[,name]),digits=2) 
  }
}
##Ordering exons
combined_final_values$exon=factor(combined_final_values$exon,levels=unique(final_values$exon),ordered=T)
combined_final_values=combined_final_values[order(combined_final_values$exon),]

## Producing final report
if (!is.null(opt$verbose)) {cat(paste("Producing final reports \n"))}

for (i in 1:length(colnames(final_values))){
  if (isTRUE(grepl("^.*mixA.*_norm_ratio$",colnames(final_values)[i],ignore.case=TRUE))) {
  samplename_report_mixA=sub("_norm_ratio","",colnames(final_values)[i])
  meancov_mixA_name=sub("_norm_ratio","_cov",colnames(final_values)[i])
  samplename_mixA=unlist(strsplit(colnames(final_values)[i],"[_]"))[5]
  samplename=unlist(strsplit(samplename_mixA,"[-]"))[1]
  #   Little fix up for case specific differences
  samplename_mixA=sub("mixA","MixA",samplename_mixA)
  samplename_mixB=sub("MixA","MixB",samplename_mixA)
  samplename=unlist(strsplit(samplename_mixA,"[-]"))[1]
  grepl_expr_mixB=paste("^.*",samplename_mixB,".*ratio",sep="")
    for (n in 1:length(colnames(final_values))) {
      if(isTRUE(grepl(grepl_expr_mixB,colnames(final_values)[n],ignore.case=TRUE))){
        meancov_mixB_name=sub("_norm_ratio","_cov",colnames(final_values)[n])
        samplename_report_mixB=sub("_norm_ratio","",colnames(final_values)[n])
        # Final report
        outfile=(paste(opt$OutputDir,"/",samplename,"_summary_report.txt",sep=''))
        sink(outfile)
        cat("Technical Report for",samplename, "\n\n")
        mean_cov_mixA=mean(final_values[,paste(meancov_mixA_name)],na.rm=TRUE)
        total_reads_mixA=sum(final_values[,paste(meancov_mixA_name)],na.rm=TRUE)
        mean_cov_mixB=mean(final_values[,paste(meancov_mixB_name)],na.rm=TRUE)
        total_reads_mixB=sum(final_values[,paste(meancov_mixB_name)],na.rm=TRUE)
        cat("Mean coverage for", samplename_report_mixA, " is: ",mean_cov_mixA,"\n")
        cat("Mean coverage for", samplename_report_mixB, " is: ",mean_cov_mixB,"\n")
        cat("Total number of aligned reads for", samplename_report_mixA, " is: ",total_reads_mixA,"\n")
        cat("Total number of aligned reads for", samplename_report_mixB, " is: ",total_reads_mixB,"\n")
        cat("Mean coverage for NTC in Mix A is:",mean_cov_NTC_mixA,"\n")
        cat("Mean coverage for NTC in Mix B is:",mean_cov_NTC_mixB,"\n\n")
        write.table(combined_final_values[,c("exon",
                              paste(samplename,"_combinedratio",sep=""),
                              paste(samplename,"_stdev_exon",sep=""),
                              paste(samplename,"_stdev_mix",sep=""))], sep="\t", row.names=F, 
              col.names=c("Exon","Mean Ratio", "StDev for exons","StDev for mixes"), quote=F)    
        sink()
        }
    }
  }
}

#Output raw coverage 
write.table(final_values[,grepl("probe|chr_pos|*_cov",colnames(final_values))],
            file=paste(opt$OutputDir,"/Raw_Coverage.txt",sep=''),sep="\t", row.names=F,quote=F)

#Plotting ratios
if (!is.null(opt$verbose)) {cat(paste("Drawing ratio plots \n"))}

for (i in 1:length(colnames(combined_final_values))) { 
  if (isTRUE(grepl("_combinedratio",colnames(combined_final_values)[i],ignore.case=TRUE))){
    name=sub("_combinedratio","",colnames(combined_final_values)[i])
    combined_final_values[,paste(name,"_stdev_combined",sep = "")]=combined_final_values[,paste(name,"_stdev_exon",sep = "")]
    combined_final_values[,paste(name,"_stdev_combined",sep = "")][is.na(combined_final_values[,paste(name,"_stdev_combined",sep = "")])]=
      combined_final_values[,paste(name,"_stdev_mix",sep = "")][is.na(combined_final_values[,paste(name,"_stdev_combined",sep = "")])]
    min_error=(combined_final_values[,colnames(combined_final_values)[i]]-combined_final_values[,paste(name,"_stdev_combined",sep = "")])
    max_error=(combined_final_values[,colnames(combined_final_values)[i]]+combined_final_values[,paste(name,"_stdev_combined",sep = "")])
    ggplot()+
      geom_point(aes(x=combined_final_values$exon,y=combined_final_values[,colnames(combined_final_values)[i]]))+
      geom_errorbar(aes(ymin=min_error,ymax=max_error,x=combined_final_values$exon),width=0)+
      geom_line(aes(x=c(1:nrow(combined_final_values)),y=0.7),colour="red")+
      geom_line(aes(x=c(1:nrow(combined_final_values)),y=1.3),colour="blue")+
      labs(x="Probes",y="Ratio (Sample to Control)",title=colnames(combined_final_values)[i])+
      theme(axis.text.x=element_text(angle=65,vjust=1,hjust=1),
            axis.title.x=element_text(size=18),
            axis.title.y=element_text(size=18),
            plot.title = element_text(size=20, face="bold", vjust=2))+
      expand_limits(x = 0, y = 0)
    ggsave(paste(opt$OutputDir,"/",colnames(combined_final_values)[i],".pdf",sep=''),height=7,width=30)
  }
}  

if (!is.null(opt$verbose)) {cat(paste(Sys.time(),"MLPAbrary Analysis finished \n"))}