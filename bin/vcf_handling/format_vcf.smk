################ CONFIGURATION ################
chroms = list(range(1,23)) + ["X", "Y"]  # No XY vars

## Merged_dataset
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
        "results/bed/merged-updated.bed"
    output:
        "results/VCF/merged_chr_{chr}.vcf"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.vcf", "", output[0])
        cmd = f"{params.plink} --real-ref-alleles --bfile {input_pattern} --chr {wildcards.chr}  --recode vcf --out {output_pattern}"
        shell(cmd)


## Single datasets
rule single_vcf:
    input:
        bed = "run_folder/update_ref/{dataset}_FORWARD.bed"
    output:
        vcf = "run_folder/single_vcf/{dataset}.vcf"
    run:
        input_pattern  = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.vcf", "", output.vcf)
        shell(f"plink1.9 --keep-allele-order --allow-no-sex --bfile {input_pattern} --recode vcf --out {output_pattern}")

