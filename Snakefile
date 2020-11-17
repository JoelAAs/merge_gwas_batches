configfile: "config.yml"
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
        "results/QC/merged.check.ref"


rule checkVCF:
    params:
        fasta_ref = config["fasta_ref"]
    input:
        "results/VCF/merged_all.vcf.gz"
    output:
          "results/QC/merged.check.af",
          "results/QC/merged.check.dup",
          "results/QC/merged.check.geno",
          "results/QC/merged.check.log",
          "results/QC/merged.check.mono",
          "results/QC/merged.check.nonSnp",
          "results/QC/merged.check.ref"
    shell:
        """
        python2.7 bin/checkVCF.py -r {params.fasta_ref} -o results/QC/merged {input}
        """


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
        cmd = f"plink1.9 --bfile {input_pattern} --chr {wildcards.chr} --recode vcf --out {output_pattern}"
        shell(cmd)

rule merge:
    """
    Merge filtered and haplotype inferred bams
    Return: Merged VCF
    """
    params:
        plink = config["plink"]
    input:
        bed = expand("run_folder/mapped/{data}_mapped.bed", data = config["datasets"]),
        fam = expand("run_folder/mapped/{data}_mapped.fam", data = config["datasets"]),
        bim = expand("run_folder/mapped/{data}_mapped.bim", data = config["datasets"])
    output:
        "results/bed/merged.bed"
    run:
        list_file = "results/merge-list.txt"
        with open(list_file, "w+") as f:
            for plink_files in zip(*[input.bed, input.bim, input.fam]):
                f.write("\t".join(plink_files) + "\n")

        output_pattern = re.sub("\\.bed", "", output[0])
        cmd = f"plink1.9 --merge-list {list_file} --make-bed --out {output_pattern}"
        shell(cmd)

# rule remove_duplicate_name:
#     """
#     Remove prefix (exm- etc) and conflicting ids within individual files
#     """
#     input:
#         mapped = "run_folder/mapped/{data}_mapped.bed"
#     output:
#         bed = "run_folder/filtered/{data}_no_prefix.bed",
#         bim = "run_folder/filtered/{data}_no_prefix.bim",
#         fam = "run_folder/filtered/{data}_no_prefix.fam"
#     shell:
#         """
#         Rscript bin/remove_duplicates_prefix.R {input.mapped} {output.bed}
#         """
#
# rule remove_duplicates_position:
#     """
#     Remove conflicting ids and positions, choosing rs-ids first other taking the majority vote
#     """
#     input:
#         mapped = get_dup_input
#     output:
#         expand("run_folder/dedup/{data}{{suffix}}dedup.bed", data = config["datasets"]),
#         expand("run_folder/dedup/{data}{{suffix}}dedup.bim", data = config["datasets"]),
#         expand("run_folder/dedup/{data}{{suffix}}dedup.fam", data = config["datasets"])
#     shell:
#         """
#         Rscript bin/remove_duplicates_position.R {wildcards.suffix} {input.mapped}
#         """

rule update_ref:
    """
    Update chromosome and position of ids given supplied dbsnp.map
    """
    input:
        get_ref_input
    log:
        "logs/{data}_update_ref.log"
    output:
        bed = "run_folder/mapped/{data}_mapped.bed",
        fam = "run_folder/mapped/{data}_mapped.fam",
        bim = "run_folder/mapped/{data}_mapped.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.bed", "", output.bed)
        strand_files = "data/update_build/InfiniumOmniExpressExome-8v1-3_A-b37-strand/InfiniumOmniExpressExome-8v1-3_A-b37.strand"
        shell(f"bin/update_build.sh {input_pattern} {strand_files} {output_pattern} {wildcards.data} 2&> {log}")

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
        shell(f"plink1.9 --file {input_pattern} --make-bed --out {output_pattern}")
