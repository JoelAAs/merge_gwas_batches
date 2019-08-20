library(snpStats)

remove_duplicate <- function(file_name, output_name) {
  file_name <- sub("\\.bed", "", file_name)
  output_name <- sub("\\.bed", "", output_name)

  plink_format <- read.plink(
    paste0(file_name, ".bed"),
    paste0(file_name, ".bim"),
    paste0(file_name, ".fam")
  )
  
  # Get duplicates
  snp_names     <- plink_format$map$snp.name 
  extra_prefix  <- sapply(snp_names[grepl(".*[_-]rs.*", snp_names)], function(x) sub("^.*-", "", x))
  remove_duplicate <- names(extra_prefix[extra_prefix %in% snp_names])
  
  # Remove duplicates with prefix and remove rest of prefixes
  plink_format$map <- plink_format$map[!plink_format$map$snp.name %in% remove_duplicate, ]
  plink_format$map$snp.name <- sapply(
    plink_format$map$snp.name, 
    function(x) sub(".*-(?=rs)", "", x,perl = T))
  row.names(plink_format$map) <- plink_format$map$snp.name
  
  # Remove duplicates from genotype matrix
  markers <- colnames(plink_format$genotypes)
  markers <- markers[!markers %in% remove_duplicate]
  plink_format$genotypes <- plink_format$genotypes[, markers]
  colnames(plink_format$genotypes) <- sapply(
    colnames(plink_format$genotypes),
    function(x) sub(".*-(?=rs)", "", x,perl = T))
  
  write.plink(
    output_name,
    snps = plink_format$genotypes,
    subject.data = plink_format$fam,
    pedigree = pedigree, id = member, father = father, mother = mother, sex = sex,
    snp.data = plink_format$map,
    chromosome = chromosome, position = position, allele.1 = allele.1, allele.2 = allele.2
    )
}

in_args     <- commandArgs(TRUE)
input_file  <- in_args[1]
output_name <- in_args[2]
remove_duplicate(input_file, output_name)
