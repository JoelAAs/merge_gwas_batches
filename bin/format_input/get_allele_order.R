# Title     : TODO
# Objective : TODO
# Created by: joel
# Created on: 2020-12-03

library(data.table)

write_allele_file <- function(manifest_file, suffix, allele_order_file) {
  flip <- function(base){
    return(
      switch(
      base,
      "A" = "T",
      "T" = "A",
      "G" = "C",
      "C" = "G",
      base
      ))
  }
  manifest_info <- fread(
    manifest_file,
    skip = "IlmnID",
    fill = TRUE)

  manifest_info <- manifest_info[!is.na(manifest_info$MapInfo)]
  manifest_info$neIlmnSource  <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 2)
  manifest_info$neDbsnpSource <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 3)

  manifest_info$Ref <- sapply(
    manifest_info$SNP,
    function(x) substring(
      strsplit(x, "/")[[1]][1],
      2
      )
    )

  manifest_info$Alt <- sapply(
    manifest_info$SNP,
    function(x) substring(
      strsplit(x, "/")[[1]][2],
      1,2
    )
  )

  if (suffix == "TOP") {
    manifest_info_flip_idx <- manifest_info$IlmnStrand == "BOT"
  } else if (suffix == "PLUS") {
    manifest_info_flip_idx <- manifest_info$RefStrand == "-"
  } else {
    exit(0)
  }
  manifest_info[manifest_info_flip_idx, "Ref"] <- sapply(
    manifest_info[manifest_info_flip_idx, ]$Ref, flip
  )

  write.table(
    manifest_info[(Chr != "" ), c('Name', "Ref")],
    col.names = FALSE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE,
    file = allele_order_file)
}


## RUN
in_args   <- commandArgs(TRUE)
manifest_path_file <- in_args[1]
suffix        <- in_args[2]
allele_order_file  <- in_args[3]

write_allele_file(manifest_path_file, suffix, allele_order_file)