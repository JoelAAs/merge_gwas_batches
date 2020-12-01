configfile: "config.yml"
import re

################ CONFIGURATION ################
chroms = list(range(1,23)) + ["X"]

wildcard_constraints:
    dataset = "[[a-zA-Z0-9-_\.]+"

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

def get_excluded_files(wildcards, suffix):
    merge_fail = "results/bed/merged-merge.missnp"
    if os.path.exists(merge_fail):
        return [f"run_folder/fliped_filtered/excluded/{data}_PLUS_updated.{suffix}" for data in config["datasets"]]
    else:
        return [f"run_folder/fliped_filtered/{data}_PLUS_updated.{suffix}" for data in config["datasets"]]


################ RULES ################

rule all:
    input:
        "results/QC/merged.check.ref",
        expand("run_folder/single_vcf/{dataset}.vcf", dataset=config["datasets"])
rule single_vcf:
    input:
        bed = "run_folder/fliped_filtered/{dataset}_PLUS_updated.bed"
    output:
        vcf = "run_folder/single_vcf/{dataset}.vcf"
    run:
        input_pattern  = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.vcf", "", output.vcf)
        shell(f"plink1.9 --bfile {input_pattern} --recode vcf --out {output_pattern}")


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
        "results/bed/merged_forward.bed"
    output:
        "results/VCF/merged_chr_{chr}.vcf"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.vcf", "", output[0])
        cmd = f"plink1.9 --bfile {input_pattern} --chr {wildcards.chr} --recode vcf --out {output_pattern}"
        shell(cmd)


rule flip_merge_to_forward:
    """
    Geno filter
    """
    input:
        bed      = "results/bed/merged.bed",
        manifest = config["final_manifest"]
    output:
        bed      = "results/bed/merged_forward.bed",
        flipfile = "run_folder/flip/merged_forward.flip",
        chrom = "run_folder/update/merged.chrom",
        pos   = "run_folder/update/merged.pos"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell("Rscript bin/update_reference.R {input.manifest} {output.chrom} 3")
        tmp = output_pattern + "_tmp"
        shell("Rscript bin/update_reference.R {input.manifest} {output.flipfile} 2")
        shell(f"plink1.9 --allow-no-sex --bfile {input_pattern} --flip {output.flipfile} --geno 0.1 --make-bed --out {tmp}")
        shell(f"plink1.9 --allow-no-sex --bfile {tmp} --extract {output.pos} --make-bed --out {output_pattern}")
        shell(f"rm {tmp}*")


rule merge:
    """
    Merge filtered bams
    Return: Merged VCF
    """
    params:
        plink = config["plink"]
    input:
        bed = lambda wildcards: get_excluded_files(wildcards, "bed"),
        fam = lambda wildcards: get_excluded_files(wildcards, "fam"),
        bim = lambda wildcards: get_excluded_files(wildcards, "bim")
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

rule exclude_misssnp:
    input:
        bed ="run_folder/fliped_filtered/{data}_PLUS_updated.bed"
    output:
        bed ="run_folder/fliped_filtered/excluded/{data}_PLUS_updated.bed",
        fam ="run_folder/fliped_filtered/excluded/{data}_PLUS_updated.fam",
        bim ="run_folder/fliped_filtered/excluded/{data}_PLUS_updated.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --bfile {input_pattern} "
              f"--exclude results/bed/merged-merge.missnp "
              f"--make-bed --out {output_pattern}")



rule update_pos_chrom:
    input:
        bed = "run_folder/fliped_filtered/{dataset}_PLUS.bed",
        chrom = "run_folder/update/{dataset}.chrom",
        pos   = "run_folder/update/{dataset}.pos"
    output:
        bed = "run_folder/fliped_filtered/{dataset}_PLUS_updated.bed",
        fam = "run_folder/fliped_filtered/{dataset}_PLUS_updated.fam",
        bim = "run_folder/fliped_filtered/{dataset}_PLUS_updated.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.bed", "", output.bed)
        tmp = input_pattern + "_tmp"
        shell(f"plink1.9 --allow-no-sex --bfile {input_pattern} --update-chr {input.chrom} --make-bed --out {tmp}_1")
        if wildcards.dataset == "TwinGene":
            shell(f"plink1.9 --allow-no-sex --bfile {tmp}_1 --update-map {input.pos} --make-bed --out {output_pattern}")
        else:
            shell(f"plink1.9 --allow-no-sex --bfile {tmp}_1 --update-map {input.pos} --make-bed --out {tmp}_2")
            shell(f"plink1.9 --allow-no-sex --bfile {tmp}_2 --extract {input.pos} --make-bed --out {output_pattern}")
        shell(f"rm {tmp}*")

rule get_pos_chrom:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        chrom = "run_folder/update/{dataset}.chrom",
        pos   = "run_folder/update/{dataset}.pos"
    shell:
        """
        Rscript bin/update_reference.R {input.manifest} {output.chrom} 3
        """

rule set_PLUS_strand:
    """
    Update chromosome and position of ids given supplied dbsnp.map
    """
    input:
        bed = "run_folder/bed/{dataset}_BOT_TOP.bed",
        flipfile = "run_folder/flip/{dataset}_PLUS.flip"
    log:
        "logs/{dataset}_update_ref.log"
    output:
        bed = "run_folder/fliped_filtered/{dataset}_PLUS.bed",
        fam = "run_folder/fliped_filtered/{dataset}_PLUS.fam",
        bim = "run_folder/fliped_filtered/{dataset}_PLUS.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --bfile {input_pattern} --flip {input.flipfile}"
              f" --not-chr 0 --set-hh-missing --make-bed --out {output_pattern} &> {log}")

rule get_flip_plus:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        flipfile = "run_folder/flip/{dataset}_PLUS.flip"
    shell:
        """
        Rscript bin/update_reference.R {input.manifest} {output.flipfile} 0
        """

rule get_flip_top_bot:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        flip_top_bot = "run_folder/flip/{dataset}_BOT_TOP.flip"
    shell:
        """
        Rscript bin/update_reference.R {input.manifest} {output.flip_top_bot} 1
        """

rule format_bed_top_bot:
    input:
        ped          = "input/{dataset}_TOP.ped",
        flip_top_bot = "run_folder/flip/{dataset}_BOT_TOP.flip"
    output:
        bed = "run_folder/bed/{dataset}_BOT_TOP.bed"
    run:
        input_pattern = re.sub("\\.ped", "", input.ped)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --file {input_pattern} --flip {input.flip_top_bot} --make-bed --out {output_pattern}")

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
