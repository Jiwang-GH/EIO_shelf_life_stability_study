# This file stores function to obtain BIP data from the Flowcentral
#Load required packages
require('RPostgreSQL')
require(tidyverse)
require(jsonlite)
require(plyr)
require(plotly)

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(colf))

## Get the path of the flowcell when providing the flowcell ID
getFCPath <- function(fcid){
  flowcentral.path <- "/ghds/flowcentral/"
  ivd.path <- "/ghds/ivd/flowcentral/"
  #return null if input is NA due to logging errors
  if (is.na(fcid)){ return(NULL)}
  flow <- list.files(path = flowcentral.path, pattern = fcid)[1]
  ivd <- list.files(path = ivd.path, pattern = fcid)[1]
  if(is.na(flow) & is.na(ivd)) {
    print(paste(fcid, "flowcell path not found!"))
    return(NULL)
  }
  if(!is.na(ivd)){
    return(paste(ivd.path, ivd, "/",sep = ""))
  }
  if(!is.na(flow)){
    return(paste(flowcentral.path, flow, "/", sep = ""))
  }
}


## QC path "/ghds/groups/reagent_qc/RQ-132_Seracare_Manual_Bip_VCV2/"
getFCPath_scvcqc <- function(fcid){
  flowcell.path <- "/ghds/groups/reagent_qc/RQ-132_Seracare_Manual_Bip_VCV2/"
  #return null if input is NA due to logging errors
  if (is.na(fcid)){ return(NULL)}
  flow <- list.files(path = flowcell.path, pattern = fcid)[1]
  if(is.na(flow)) {
    print(paste(fcid, "flowcell path not found!"))
    return(NULL)
  }
  if(!is.na(flow)){
    return(paste(flowcell.path, flow, "/", sep = ""))
  }
}

## QC path for old flowcell "/ghds/groups/bioinformatics_research/flowcentral/"
getFCPath_oldflowcell <- function(fcid){
  oldflowcell.path <- "/ghds/groups/bioinformatics_research/flowcentral/"
  #return null if input is NA due to logging errors
  if (is.na(fcid)){ return(NULL)}
  flow <- list.files(path = oldflowcell.path, pattern = fcid)[1]
  if(is.na(flow)) {
    print(paste(fcid, "flowcell path not found!"))
    return(NULL)
  }
  if(!is.na(flow)){
    return(paste(oldflowcell.path, flow, "/", sep = ""))
  }
}



## Pull flowcell sequencing metrics
# interop_db.hdr.tsv
pull_seq_info <- function(fcid_list){
  ser_final <- list()
  for(x in 1:length(fcid_list)){
    # x <- 1
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    seq_file <- list.files(fcid, pattern = glob2rx("interop_db.hdr.tsv"), full.names = TRUE)
    seq_list <- sapply(seq_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    seq_list <- lapply(seq_list, transform)
    seq <- bind_rows(seq_list, .id = 'file')
    seq_final <- rbind(seq_final, seq)
  }
  return(seq_final)
}

##Pull data about GCIQR from flowcentral based on sample ID 
pull_gciqr_info <- function(sampleid_list){
  gciqr_final <- list()
  Fpath1 <- "/ghds/ivd/flowcentral/"
  Fpath2 <- "/ghds/flowcentral/"
  for(x in 1:length(sampleid_list)){
    sample <- sampleid_list[x]
    fileindex <- paste(sample, "*.ghcnv_qc.hdr.tsv", sep = "")
    ##Two line of file used for the GC IQR analysis
    gciqr_file1 <- Sys.glob(file.path(Fpath1, "[0-9]*[A-Z]", fileindex))
    gciqr_file2 <- Sys.glob(file.path(Fpath2, "[0-9]*[A-Z]", fileindex))
    ## New analysis should use this one,the matching pattern is better
    # gciqr_file1 <- Sys.glob(file.path(Fpath1,  "[0-9]*[A-Z][^.]", fileindex))
    # gciqr_file2 <- Sys.glob(file.path(Fpath2,  "[0-9]*[A-Z][^.]", fileindex))
   
    gciqr_file <- c(gciqr_file1, gciqr_file2)
    gciqr_list <- sapply(gciqr_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    gciqr_list <- lapply(gciqr_list, transform, run_sample_id = as.character(run_sample_id))
    gciqr <- bind_rows(gciqr_list, .id = 'file')
    gciqr_final <- rbind(gciqr_final, gciqr)
  }
  return(gciqr_final)
}

##Pull data about GCIQR from flowcentral based on sample ID 
pull_gciqr_info1 <- function(sampleid_list){
  gciqr_final <- list()
  Fpath1 <- "/ghds/ivd/flowcentral/"
  Fpath2 <- "/ghds/flowcentral/"
  for(x in 1:length(sampleid_list)){
    sample <- sampleid_list[x]
    fileindex <- paste(sample, "*.ghcnv_qc.hdr.tsv", sep = "")
    gciqr_file1 <- Sys.glob(file.path(Fpath1,  "[0-9]*[A-Z][^.]", fileindex))
    gciqr_file2 <- Sys.glob(file.path(Fpath2,  "[0-9]*[A-Z][^.]", fileindex))
    
    gciqr_file <- c(gciqr_file1, gciqr_file2)
    gciqr_list <- sapply(gciqr_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    gciqr_list <- lapply(gciqr_list, transform, run_sample_id = as.character(run_sample_id))
    gciqr <- bind_rows(gciqr_list, .id = 'file')
    gciqr_final <- rbind(gciqr_final, gciqr)
  }
  return(gciqr_final)
}

##Pull the insertsize data from BIP output 
pull_isize_info <- function(fcid_list){
  isize_final  <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    isize_file <- list.files(fcid, pattern = glob2rx("*.short.bwamem.cram.isize.txt"), full.names = TRUE)
    isize_list <- sapply(isize_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    isize_list <- lapply(isize_list, transform)
    isize <- bind_rows(isize_list, .id = 'file')
    isize_final <- rbind(isize_final, isize)
  }
  return(isize_final)
}

# Pull the coverage informaiton from BIP output 
pull_coverage_info <- function(fcid_list){
  coverage_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    coverage_file <- list.files(fcid, pattern = glob2rx("*.coverage.hdr.tsv"), full.names = TRUE)
    coverage_list <- sapply(coverage_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    coverage_list <- lapply(coverage_list, transform, run_sample_id = as.character(run_sample_id) )
    coverage <- rbind.fill(coverage_list,  .id = 'file')
    coverage_final <- rbind(coverage_final, coverage)
  }
  return(coverage_final)
}


# Pull the coverage informaiton from BIP output for seracare path
pull_coverage_info_scvcqc <- function(fcid_list){
  coverage_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath_scvcqc(as.character(fcid_list[[x]]))
    coverage_file <- list.files(fcid, pattern = glob2rx("*.coverage.hdr.tsv"), full.names = TRUE)
    coverage_list <- sapply(coverage_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    coverage_list <- lapply(coverage_list, transform, run_sample_id = as.character(run_sample_id) )
    coverage <- rbind.fill(coverage_list,  .id = 'file')
    coverage_final <- rbind(coverage_final, coverage)
  }
  return(coverage_final)
}

# Pull the autoqc informaiton from BIP output 
pull_autoqc_info_scvcqc <- function(fcid_list){
  autoqc_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath_scvcqc(as.character(fcid_list[[x]]))
    autoqc_file <- list.files(fcid, pattern = glob2rx("autoqc_sample_qc.hdr.tsv"), full.names = TRUE)
    autoqc_list <- sapply(autoqc_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    autoqc_list <- lapply(autoqc_list, transform, run_sample_id = as.character(run_sample_id) )
    autoqc <- rbind.fill(autoqc_list,  .id = 'file')
    autoqc_final <- rbind(autoqc_final, autoqc)
  }
  return(autoqc_final)
}


# Pull the autoqc informaiton from BIP output 
pull_autoqc_info <- function(fcid_list){
  autoqc_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    autoqc_file <- list.files(fcid, pattern = glob2rx("autoqc_sample_qc.hdr.tsv"), full.names = TRUE)
    autoqc_list <- sapply(autoqc_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    autoqc_list <- lapply(autoqc_list, transform, run_sample_id = as.character(run_sample_id) )
    autoqc <- rbind.fill(autoqc_list,  .id = 'file')
    autoqc_final <- rbind(autoqc_final, autoqc)
  }
  return(autoqc_final)
}

## Pull ontarget db data from BIP output
# .on_target_db.hdr.tsv
pull_otdb_info <- function(fcid_list){
  otdb_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    otdb_file <- list.files(fcid, pattern = glob2rx("*.on_target_db.hdr.tsv"), full.names = TRUE)
    otdb_list <- sapply(otdb_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    otdb_list <- lapply(otdb_list, transform, run_sample_id = as.character(run_sample_id))
    otdb <- bind_rows(otdb_list, .id = 'file')
    otdb_final <- rbind(otdb_final, otdb)
  }
  return(otdb_final)
}

# Pull probe infor from BIP output
pull_probe_info <- function(fcid_list){
  probe_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    probe_file <- list.files(fcid, pattern = glob2rx("*.ghcnv_probe.hdr.tsv"), full.names = TRUE)
    probe_list <- sapply(probe_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    probe_list <- lapply(probe_list, transform)
    probe <- bind_rows(probe_list, .id = 'file')
    probe_final <- rbind(probe_final, probe  )
  }
  return( probe_final)
}


## Pull the HS/BB ratio data from BIP output from a QC path "/ghds/groups/reagent_qc/RQ-132_Seracare_Manual_Bip_VCV2/"

pull_hsbb_info_qc <- function(fcid_list){
  hss_bb_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath_scvcqc(as.character(fcid_list[[x]]))
    hsbb_file <- list.files(fcid, pattern = glob2rx("*.ghcnv_qc.hdr.tsv"), full.names = TRUE)
    hsbb_list <- sapply(hsbb_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    hsbb_list <- lapply(hsbb_list, transform, run_sample_id = as.character(run_sample_id))
    hss_bb <- bind_rows(hsbb_list, .id = 'file')
    hss_bb_final <- rbind(hss_bb_final, hss_bb)
  }
  return(hss_bb_final)
}

## Pull the HS/BB ratio data from BIP output
pull_hsbb_info <- function(fcid_list){
  hss_bb_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    hsbb_file <- list.files(fcid, pattern = glob2rx("*.ghcnv_qc.hdr.tsv"), full.names = TRUE)
    hsbb_list <- sapply(hsbb_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    hsbb_list <- lapply(hsbb_list, transform, run_sample_id = as.character(run_sample_id))
    hss_bb <- bind_rows(hsbb_list, .id = 'file')
    hss_bb_final <- rbind(hss_bb_final, hss_bb)
  }
  return(hss_bb_final)
}

## Pull MSI sites from BIP output
pull_msi_sites = function(fcid_list){
  msi_final<- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath_scvcqc(as.character(fcid_list[x]))
    msi_file <- list.files(fcid, pattern = glob2rx("*.msi_parameters.hdr.tsv"), full.names = TRUE)
    msi_list <- sapply(msi_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    msi_list <- lapply(msi_list, transform, run_sample_id = as.character(run_sample_id) )
    msi <- rbind.fill(msi_list,  .id = 'file')
    msi_final <- rbind(msi_final, msi)
  }
  return(msi_final)
}

## pull SNV variant from BIP output
pull_variant = function(fcid_list){
variant_final <- list()
for(x in 1:length(fcid_list)){
  fcid <- getFCPath_scvcqc(as.character(fcid_list[x]))
  variant_file <- list.files(fcid, pattern = glob2rx("*gh_variant.hdr.tsv"), full.names = TRUE)
  variant_list <- sapply(variant_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
  variant_list <- lapply(variant_list, transform, run_sample_id = as.character(run_sample_id))
  variant_list <- lapply(variant_list, transform, chrom = as.character(chrom))
  variant <- bind_rows(variant_list, .id = 'file')
  variant_final <- rbind( variant_final, variant)
}
return(variant_final)
}



## pull indel variant from BIP output
pull_indel_call = function(fcid_list){
  indel_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath_scvcqc(as.character(fcid_list[[x]]))
    indel_file <- list.files(fcid, pattern = glob2rx("*gh_indel.hdr.tsv"), full.names = TRUE)
    indel_list <- sapply(indel_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    indel_list <- lapply(indel_list, transform, run_sample_id = as.factor(run_sample_id))
    indel_list <- lapply(indel_list, transform, run_sample_id = as.character(run_sample_id))
    indel <- bind_rows(indel_list, .id = 'file')
    indel_final <- rbind(indel_final, indel)
  }
  return(indel_final)
}


## Pull CNV call from BIP output
pull_cnv_call = function(fcid_list){
  cnv_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath(as.character(fcid_list[[x]]))
    cnv_file <- list.files(fcid, pattern = glob2rx("*cnv_call.hdr.tsv"), full.names = TRUE)
    cnv_list <- sapply(cnv_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    cnv_list <- lapply(cnv_list, transform, run_sample_id = as.factor(run_sample_id))
    cnv_list <- lapply(cnv_list, transform, chrom = as.character(chrom))
    cnv <- bind_rows(cnv_list, .id = 'file')
    cnv_final <- rbind(cnv_final, cnv)
  }
  return(cnv_final)
}

## Pull fusion call from BIP output
pull_fusion_call = function(fcid_list){
  fusion_final <- list()
  for(x in 1:length(fcid_list)){
    fcid <- getFCPath_scvcqc(as.character(fcid_list[[x]]))
    fusion_file <- list.files(fcid, pattern = glob2rx("*.fusion_call.hdr.tsv"), full.names = TRUE)
    fusion_list <- sapply(fusion_file, read.delim, simplify = FALSE, USE.NAMES = TRUE)
    fusion_list <- lapply(fusion_list, transform, run_sample_id = as.character(run_sample_id))
    fusion <- bind_rows(fusion_list, .id = 'file')
    fusion_final <- rbind(fusion_final, fusion)
  }
  return(fusion_final)
}

## the theme_Publication
theme_Publication <- function(base_size=20, base_family="helvetica") {
# theme_Publication <- function(base_size, base_family) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

## Pull input file to get the cfDNA lot # and Input for each of the sample
pull_vc_input_files <- function(){
  Fpath <- "/ghds/groups/VC_characterization/INPUT_files/"
  input_file <-  list.files(Fpath, pattern= "*.csv", include.dirs = TRUE, full.names = TRUE)
  inputfile_data <- rbindlist(lapply(input_file,fread))
  return(inputfile_data)
}



