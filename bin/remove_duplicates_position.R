library(snpStats)
library(dplyr)
library(tidyverse)
library(parallel)

# Wrapper for read.plink from plink collection path
read_plink_wrap <- function(file_name) {
  file_name <- sub("\\.bed", "", file_name)
  
  plink_format <- read.plink(
    paste0(file_name, ".bed"),
    paste0(file_name, ".bim"),
    paste0(file_name, ".fam")
  )
  return(plink_format)
}
 

# Return Dataframe with coulmns "snp.name" and "position" from all input bed files
read_snp_pos_data <- function(plink_list){
  columns <- c("snp.name", "position")
  df_pos_snp_list <- lapply(plink_list, function(x) x$map[, columns])
  for (df in df_pos_snp_list){
    row.names(df) <- NULL
  }
  return(bind_rows(df_pos_snp_list))
}

mark_duplicates <- function(df_pos_snp){
  print("marking duplicates")
  cl <- makeCluster(detectCores() - 1)
  df_pos_snp              <- df_pos_snp[!is.na(df_pos_snp$position), ]
  df_pos_snp$concat       <- parApply(cl, df_pos_snp, 1, function(x) paste0(trimws(x[1]) ,trimws(x[2])))
  stopCluster(cl)
  df_pos_snp$snp_count    <- unname(table(df_pos_snp$snp.name)[df_pos_snp$snp.name])
  df_pos_snp$pos_count    <- unname(table(df_pos_snp$position)[as.character(df_pos_snp$position)])
  df_pos_snp$concat_count <- unname(table(df_pos_snp$concat)[df_pos_snp$concat])
  df_pos_snp %>%  
    filter(snp_count !=  pos_count | pos_count != concat_count) %>% 
    return()
}


remove_duplicates <- function(multi_pos_df, snp_name=T) {
  print("Removing duplicates")
  if (!snp_name) {
    multi_name <- multi_pos_df %>%
      filter(pos_count > concat_count & snp_count == concat_count)
  } else{
    multi_name <- multi_pos_df %>%
      filter(snp_count > concat_count & pos_count == concat_count)
  }
  
  id_to_keep <- function(current_pos, multi_name) {
    current_lines = multi_name[multi_name$position == current_pos, ]
    snp_table = table(current_lines$snp.name)
    rs_table = snp_table[grepl("^rs.*", names(snp_table))]
    if (length(rs_table) >= 1) {
      keep_id = names(which.max(rs_table))
    } else {
      keep_id = names(which.max(snp_table))
    }
    return(keep_id)
  }
  cl <- makeCluster(detectCores() - 1)
  print(head(multi_name))
  clusterExport(cl, varlist =  list("multi_name", "id_to_keep"), envir = environment())
  unique_pos = unique(multi_name$position)
  rsid_keep <- parSapply(cl = cl, unique_pos, FUN =  function(x) id_to_keep(x, multi_name))
  stopCluster(cl)
  
  multi_pos_df %>%  filter(!snp.name %in% rsid_keep) %>% 
    {unique(.$snp.name)} %>% 
    return()
}

write_to_output <- function(rsid_remove, plinks, suffix){
  print("Writing to output")
  names_data = names(plinks)
  print(names_data)

  for (i in 1:length(plinks)) {
    plink_format = plinks[[i]]
    output_name = paste0("run_folder/dedup/", names_data[[i]], suffix ,"dedup")
    
    na_pos_id  <- plink_format$map[is.na(plink_format$map$position), "snp.name"]
    all_remove <- append(rsid_remove, na_pos_id)
    
    plink_format$map <- plink_format$map[!plink_format$map$snp.name %in% all_remove, ]
    
    markers <- colnames(plink_format$genotypes)
    markers <- markers[!markers %in% all_remove]
    plink_format$genotypes <- plink_format$genotypes[, markers]
    
    write.plink(
      output_name,
      snps = plink_format$genotypes,
      subject.data = plink_format$fam,
      pedigree = pedigree, id = member, father = father, mother = mother, sex = sex,
      snp.data = plink_format$map,
      chromosome = chromosome, position = position, allele.1 = allele.1, allele.2 = allele.2
    )
  } 
}

#### RUN
in_args     <- commandArgs(TRUE)
input_files <- in_args[2:length(in_args)]
suffix <- in_args[1]

input_names <- unname(sapply(input_files, function(x) sub("(_no_prefix_mapped|_flipped)\\.bed", "", tail(strsplit(x, "/")[[1]], 1))))
plinks <- lapply(input_files, read_plink_wrap)
names(plinks) <- input_names

all_df_pos_snp <- read_snp_pos_data(plinks)

remove_pos <- all_df_pos_snp %>% 
  mark_duplicates() %>% 
  remove_duplicates(snp_name=F)

remove_snp_pos <- all_df_pos_snp %>%
  filter(!snp.name %in% remove_pos) %>% 
  mark_duplicates() %>% 
  remove_duplicates(snp_name=T)

append(remove_pos, remove_snp_pos) %>% 
  write_to_output(plinks = plinks, suffix = suffix)

