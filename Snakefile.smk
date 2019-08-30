configfile: "config.yml"
import os
import re

################ CONFIGURATION ################
chroms = list(range(1,23)) + ["X"]

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
        "results/VCF/merged_all.vcf.gz"

rule merge_chr:
    input:
        expand("results/VCF/merged_chr_{chr}.vcf.gz", chr=chroms)
    output:
        "results/VCF/merged_all.vcf.gz"
    shell:
        """
        vcf-concat {input} | gzip -c > {output}
        """

rule sort_and_zip:
    input:
        "results/VCF/merged_chr_{chr}.vcf"
    output:
        "results/VCF/merged_chr_{chr}.vcf.gz"
    shell:
        """
        vcf-sort {input} | bgzip -c > {output}
        """

rule recode_vcf_per_chrom:
    params:
        plink = config["plink"]
    input:
        "results/bed/merged.bed"
    output:
        "results/VCF/merged_chr_{chr}.vcf"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.vcf", "", output[0])
        cmd = f"{params.plink} --bfile {input_pattern} --mind 0.1 --chr {wildcards.chr} --recode vcf --out {output_pattern}"
        print(cmd)
        shell(cmd)


rule merge:
    """
    Merge filtered and haplotype inferred bams
    Return: Merged VCF
    """
    params:
        plink = config["plink"]
    input:
        bed = expand("run_folder/dedup/{data}_dedup.bed", data = config["datasets"]),
        fam = expand("run_folder/dedup/{data}_dedup.fam", data = config["datasets"]),
        bim = expand("run_folder/dedup/{data}_dedup.bim", data = config["datasets"])
    output:
        "results/bed/merged.bed"
    run:
        list_file = "results/merge-list.txt"
        with open(list_file, "w+") as f:
            for plink_files in zip(*[input.bed[1:], input.bim[1:], input.fam[1:]]):
                f.write("\t".join(plink_files) + "\n")

        first_input = re.sub("\\.bed", "", input.bed[0])
        output_pattern = re.sub("\\.bed", "", output[0])
        cmd = f"{params.plink} --bfile {first_input} --merge-list {list_file} --make-bed --out {output_pattern}"
        shell(cmd)


rule remove_duplicate_name:
    """
    Remove prefix (exm- etc) and conflicting ids within individual files
    """
    input:
        mapped = "run_folder/mapped/{data}_mapped.bed"
    output:
        bed = "run_folder/filtered/{data}_no_prefix.bed",
        bim = "run_folder/filtered/{data}_no_prefix.bim",
        fam = "run_folder/filtered/{data}_no_prefix.fam"
    shell:
        """
        Rscript bin/remove_duplicates_prefix.R {input.mapped} {output.bed}
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
