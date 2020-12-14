## Input function
def get_excluded_files(wildcards, suffix):
    merge_fail = "results/bed/merged-merge.missnp"
    if os.path.exists(merge_fail):
        return [f"run_folder/update_ref/excluded/{data}_FORWARD.{suffix}" for data in config["datasets"]]
    else:
        return [f"run_folder/update_ref/{data}_FORWARD-updated.{suffix}" for data in config["datasets"]]


## Rules
rule run_plink_update:
    input:
        "results/bed/Run-plink.sh"
    output:
        "results/bed/merged-updated.bed"
    shell:
        """
        plink1.9 --bfile results/bed/merged --exclude results/bed/Exclude-merged-1000G.txt --make-bed --out results/bed/TEMP1
        plink1.9 --bfile results/bed/TEMP1 --update-map results/bed/Chromosome-merged-1000G.txt --update-chr --make-bed --out results/bed/TEMP2
        plink1.9 --bfile results/bed/TEMP2 --update-map results/bed/Position-merged-1000G.txt --make-bed --out results/bed/TEMP3
        plink1.9 --bfile results/bed/TEMP3 --flip results/bed/Strand-Flip-merged-1000G.txt --make-bed --out results/bed/TEMP4
        plink1.9 --bfile results/bed/TEMP4 --a2-allele results/bed/Force-Allele1-merged-1000G.txt --make-bed --out results/bed/merged-updated
        rm results/bed/TEMP*
        """

rule update_1000g:
    params:
        hrc_check = "/home/joel/Documents/software/HRC-1000G-check-bim/HRC-1000G-check-bim.pl",
        ref_panel_1000g = "data/1000GP_Phase3_combined.legend"
    input:
        bim  = "results/bed/merged.bim",
        freq = "results/bed/freq/merged.frq"
    output:
        "results/bed/Run-plink.sh"
    shell:
        """
        perl {params.hrc_check} -b {input.bim} -f {input.freq} -r {params.ref_panel_1000g} -g -p EUR -s Run-plink.sh
        """


rule get_frequency:
    input:
        bed = "results/bed/merged.bed"
    output:
        freq = "results/bed/freq/merged.frq"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.frq", "", output.freq)
        shell("plink1.9 --bfile {input_pattern} --real-ref-alleles --freq --out {output_pattern}")

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
        "results/bed/merged.bed",
        "results/bed/merged.bim"
    run:
        list_file = "results/merge-list.txt"
        with open(list_file, "w+") as f:
            for plink_files in zip(*[input.bed, input.bim, input.fam]):
                f.write("\t".join(plink_files) + "\n")

        output_pattern = re.sub("\\.bed", "", output[0])
        cmd = f"plink1.9 --keep-allele-order --allow-no-sex --merge-list {list_file} --geno --maf --make-bed --out {output_pattern}"
        shell(cmd)


rule exclude_misssnp:
    input:
        bed ="run_folder/update_ref/{data}_FORWARD.bed",
        misssnp = "results/bed/merged-merge.missnp"
    output:
        bed ="run_folder/update_ref/excluded/{data}_FORWARD.bed",
        fam ="run_folder/update_ref/excluded/{data}_FORWARD.fam",
        bim ="run_folder/update_ref/excluded/{data}_FORWARD.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --keep-allele-order --allow-no-sex --bfile {input_pattern} "
              f"--exclude {input.misssnp} "
              f"--make-bed --out {output_pattern}")
