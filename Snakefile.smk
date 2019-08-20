configfile: "config.yml"
import os
import re

################ CONFIGURATION ################
chroms = list(range(1,23)) + ["x", "y"]

################ FUNCTIONS ################
def get_ref_input(wildcards):
    """
    Choose input given wildcard for rule "update_ref"
    """
    input_pattern = wildcards["data"]
    if "_no_prefix" in input_pattern:
        return f"run_folder/filtered/{input_pattern}.bed"
    else:
        return f"run_folder/bed/{input_pattern}.bed"

def get_dup_input(wildcards):
    """
    Choose input given wildcard for rule "remove_duplicates_position"
    """
    data_pattern = wildcards["suffix"]
    if "_flipped" in data_pattern:
        input_pattern = f"run_folder/realign/{{data}}_full_flipped.bed"
        return [input_pattern.format(data = dataset) for dataset in config["datasets"]]
    else:
        input_pattern = f"run_folder/mapped/{{data}}_no_prefix_mapped.bed"
        return [input_pattern.format(data = dataset) for dataset in config["datasets"]]


################ RULES ################

rule all:
    input:
        #"results/merged.bed"
        expand("run_folder/realign/{data}_chrom_{chr}.phased.haps", data = config["datasets"], chr=chroms)

rule merge:
    """
    Merge filtered and haplotype inferred bams
    Return: Merged VCF
    """
    params:
        plink = config["plink"]
    input:
        bed = expand("run_folder/dedup/{data}_full_flipped_dedup.bed", data = config["datasets"]),
        fam = expand("run_folder/dedup/{data}_full_flipped_dedup.fam", data = config["datasets"]),
        bim = expand("run_folder/dedup/{data}_full_flipped_dedup.bim", data = config["datasets"])
    output:
        "results/merged.bed"
    shell:
        """
        echo "hej"
        """

rule flip_join:
    params:
        plink = config["plink"]
    input:
        expand("run_folder/realign/{{data}}_chrom_{chr}_flipped.bed",  chr=chroms)
    output:
        "run_folder/realign/{data}_full_flipped.bed",
    shell:
        """
        """

rule flip:
    """
    Flip variants according to haplotype definitions
    """
    params:
        plink = config["plink"]
    input:
        bed = "run_folder/splits/{data}_chrom_{chr}.bed",
        bim = "run_folder/splits/{data}_chrom_{chr}.bim",
        fam = "run_folder/splits/{data}_chrom_{chr}.fam",
        haps = "run_folder/realign/{data}_chrom_{chr}.phased.haps",
        sample = "run_folder/realign/{data}_chrom_{chr}.phased.sample"
    output:
        bed = "run_folder/realign/{data}_chrom_{chr}_flipped.bed",
        fam = "run_folder/realign/{data}_chrom_{chr}_flipped.fam",
        bim = "run_folder/realign/{data}_chrom_{chr}_flipped.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"{params.plink} --bfile {input_pattern} --flip  --make-bed --out {output_pattern}") #TODO: Psuedo 


rule realign_shapeit:
    """
    Get infered haplotype definitions supplied by dbsnp
    """
    params:
        dbsnp_map = config["dbsnp_map"],
        shapeit = config["shapeit"]
    input:
        bed = "run_folder/splits/{data}_chrom_{chr}.bed",
        bim = "run_folder/splits/{data}_chrom_{chr}.bim",
        fam = "run_folder/splits/{data}_chrom_{chr}.fam"
    output:
        haps = "run_folder/realign/{data}_chrom_{chr}.phased.haps",
        sample = "run_folder/realign/{data}_chrom_{chr}.phased.sample"
    shell:
        """
        {params.shapeit} -B {input.bed} {input.bim} {input.fam} \
            -M {params.dbsnp_map} \
            --output-max {output.haps} {output.sample} \
            --thread 2
        """

rule split_per_chromosome:
    """
    """
    params:
        plink     = config["plink"]
    input:
        bam = "run_folder/dedup/{data}_dedup.bed",
        bim = "run_folder/dedup/{data}_dedup.fam",
        fam = "run_folder/dedup/{data}_dedup.bim"
    output:
        bam = expand("run_folder/splits/{{data}}_chrom_{chr}.bed", chr=chroms),
        bim = expand("run_folder/splits/{{data}}_chrom_{chr}.bim", chr=chroms),
        fam = expand("run_folder/splits/{{data}}_chrom_{chr}.fam", chr=chroms)
    run:
        input_pattern = re.sub("\\.bed", "", input.bam)
        for chr in chroms:
            output_pattern = f"run_folder/splits/{wildcards.data}_chrom_{chr}"
            shell(f"{params.plink} --bfile {input_pattern} --chr {chr} --recode --make-bed --out {output_pattern}")

rule remove_duplicate_name:
    """
    Remove prefix (exm- etc) and conflicting ids within individual files
    """
    input:
        mapped = "run_folder/mapped/{data}_mapped.bed"
    output:
        filtered = "run_folder/filtered/{data}_no_prefix.bed"
    shell:
        """
        Rscript bin/remove_duplicates_prefix.R {input.mapped} {output.filtered}
        """

rule remove_duplicates_position:
    """
    Remove conflicting ids and positions, choosing rs-ids first other taking the majority vote
    """
    input:
        mapped = get_dup_input
    output:
        expand("run_folder/dedup/{data}{{suffix}}dedup.bed", data = config["datasets"]),
        expand("run_folder/dedup/{data}{{suffix}}dedup.bim", data = config["datasets"]),
        expand("run_folder/dedup/{data}{{suffix}}dedup.fam", data = config["datasets"])
    shell:
        """
        Rscript bin/remove_duplicates_position.R {wildcards.suffix} {input.mapped}
        """

rule update_ref:
    """
    Update chromosome and position of ids given supplied dbsnp.map
    """
    params:
        dbsnp_map = config["dbsnp_map"],
        dbsnp_chr = config["dbsnp_chr"],
        plink     = config["plink"]
    input:
        get_ref_input
    output:
        bed = "run_folder/mapped/{data}_mapped.bed"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"{params.plink} --bfile {input_pattern} --update-chr {params.dbsnp_chr} --make-bed --out {output_pattern}_chr")
        shell(f"{params.plink} --bfile {output_pattern}_chr --update-map {params.dbsnp_map} --make-bed --out {output_pattern}")
        shell("rm {output_pattern}_chr.*")

rule format_bed:
    """
    Format pad to bed
    """
    params:
        plink     = config["plink"]
    input:
        ped = "input/{dataset}.ped"
    output:
        bed = "run_folder/bed/{dataset}.bed"
    run:
        input_pattern = re.sub("\\.ped", "", input.ped)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"{params.plink} --file {input_pattern} --make-bed --out {output_pattern}")
