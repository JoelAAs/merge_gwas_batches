# Title     : Get_fliplist
# Objective : Get variants to flip
# Created by: Joel
# Created on: 2020-11-17

library(data.table)

flip_to_forward_from_plus <- function(file_path, output) {
  manifest_info <- fread(
    file_path,
    skip = "IlmnID",
    fill = TRUE)
  manifest_info$neIlmnSource  <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 2)
  manifest_info$neDbsnpSource <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 3)
  manifest_info[,
    flip := ifelse((RefStrand == '-' & neDbsnpSource == 'F') | (RefStrand == '+' & neDbsnpSource == 'R') , 1,0)
  ]
  manifest_info[flip==1,c('IlmnID','IlmnStrand','neIlmnSource','neDbsnpSource','SNP','RefStrand','flip')]
  print(paste0("Found: ", sum(manifest_info$flip), " variants to flip." ))
  write.table(
    manifest_info[flip==1,'Name'],
    col.names = FALSE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE,
    file = output)
}

flip_to_plus <- function(file_path, output) {
  manifest_info <- fread(
    file_path,
    skip = "IlmnID",
    fill = TRUE)
  manifest_info$neIlmnSource  <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 2)
  manifest_info$neDbsnpSource <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 3)
  manifest_info[,
    flip := ifelse(RefStrand == '-' , 1,0)
  ]
  manifest_info[flip==1,c('IlmnID','IlmnStrand','neIlmnSource','neDbsnpSource','SNP','RefStrand','flip')]
  print(paste0("Found: ", sum(manifest_info$flip), " variants to flip." ))
  write.table(
    manifest_info[flip==1,'Name'],
    col.names = FALSE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE,
    file = output)
}


flip_to_BOT_TOP <- function(file_path, output) {
  manifest_info <- fread(
    file_path,
    skip = "IlmnID",
    fill = TRUE)
  manifest_info$neIlmnSource  <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 2)
  manifest_info$neDbsnpSource <- sapply(strsplit(manifest_info$IlmnID,'_'), "[", 3)
  manifest_info[,
    flip := ifelse(IlmnStrand == "BOT", 1,0)
  ]

  write.table(
    manifest_info[flip==1,'Name'],
    col.names = FALSE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE,
    file = output)

}

get_pos_chrom <- function(file_path, output) {
  manifest_info <- fread(
    file_path,
    skip = "IlmnID",
    fill = TRUE)

  manifest_info <- manifest_info[(Chr != "" ), ]
  # manifest_info_idx <- manifest_info[, "Chr"] %in% c("X", "XY")
  # manifest_info[manifest_info_idx, "Chr"] <- "23"
  write.table(
    manifest_info[(Chr != "" ), c('Name', "Chr")],
    col.names = FALSE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE,
    file = output)

  output_pos <- gsub(".chrom", ".pos", output)
  write.table(
    manifest_info[!is.na(MapInfo), c('Name', "MapInfo")],
    col.names = FALSE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE,
    file = output_pos)

}


## RUN
in_args   <- commandArgs(TRUE)
file_path <- in_args[1]
output    <- in_args[2]
all_TOP   <- in_args[3]

if (all_TOP == "0") {
  print("Flip list for forward strand")
  flip_to_plus(file_path, output)
} else if (all_TOP == "1") {
  flip_to_BOT_TOP(file_path, output)
} else if (all_TOP == "2"){
  flip_to_forward_from_plus(file_path, output)
} else {
  get_pos_chrom(file_path, output)
}