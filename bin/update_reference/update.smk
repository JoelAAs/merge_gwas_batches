rule run_plink_update_single:
    input:
        "run_folder/update_ref/{dataset}.sh"
    output:
        bed="run_folder/update_ref/{dataset}_FORWARD-updated.bed",
        bim="run_folder/update_ref/{dataset}_FORWARD-updated.bim",
        fam="run_folder/update_ref/{dataset}_FORWARD-updated.fam"
    shell:
        """
        plink1.9 --bfile run_folder/update_ref/{wildcards.dataset}_FORWARD  \
            --exclude run_folder/update_ref/Exclude-{wildcards.dataset}_FORWARD-1000G.txt \
            --make-bed --out run_folder/update_ref/TEMP1_{wildcards.dataset}
        plink1.9 --bfile run_folder/update_ref/TEMP1_{wildcards.dataset} \
            --update-chr run_folder/update_ref/Chromosome-{wildcards.dataset}_FORWARD-1000G.txt \
            --make-bed --out run_folder/update_ref/TEMP2_{wildcards.dataset}
        plink1.9 --bfile run_folder/update_ref/TEMP2_{wildcards.dataset} \
            --update-map run_folder/update_ref/Position-{wildcards.dataset}_FORWARD-1000G.txt --make-bed \
            --out run_folder/update_ref/TEMP3_{wildcards.dataset}
        plink1.9 --bfile run_folder/update_ref/TEMP3_{wildcards.dataset} \
            --flip run_folder/update_ref/Strand-Flip-{wildcards.dataset}_FORWARD-1000G.txt \
            --make-bed --out run_folder/update_ref/TEMP4_{wildcards.dataset}
        plink1.9 --bfile run_folder/update_ref/TEMP4_{wildcards.dataset} \
            --a2-allele run_folder/update_ref/Force-Allele1-{wildcards.dataset}_FORWARD-1000G.txt \
            --make-bed --out run_folder/update_ref/{wildcards.dataset}_FORWARD-updated
        rm run_folder/update_ref/TEMP*_{wildcards.dataset}*
        """

rule update_1000g_single:
    params:
        hrc_check = "/home/joel/Documents/software/HRC-1000G-check-bim/HRC-1000G-check-bim.pl",
        ref_panel_1000g = "data/1000GP_Phase3_combined.legend"
    input:
        bim  = "run_folder/update_ref/{dataset}_FORWARD.bim",
        freq = "run_folder/freq/{dataset}_FORWARD.frq"
    output:
        "run_folder/update_ref/{dataset}.sh"
    shell:
        """
        perl {params.hrc_check} -b {input.bim} -f {input.freq} -r {params.ref_panel_1000g} -g -p EUR -s {wildcards.dataset}.sh
        """


rule get_frequency_single:
    input:
        bed = "run_folder/update_ref/{dataset}_FORWARD.bed"
    output:
        freq = "run_folder/freq/{dataset}_FORWARD.frq"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.frq", "", output.freq)
        shell("plink1.9 --bfile {input_pattern} --real-ref-alleles --freq --out {output_pattern}")


rule update_pos_chrom:
    input:
        bed = "run_folder/bed/{dataset}_FORWARD.bed",
        chrom = "run_folder/update_files/{dataset}.chrom",
        pos   = "run_folder/update_files/{dataset}.pos"
    output:
        bed ="run_folder/update_ref/{dataset}_FORWARD.bed",
        fam ="run_folder/update_ref/{dataset}_FORWARD.fam",
        bim ="run_folder/update_ref/{dataset}_FORWARD.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.bed", "", output.bed)
        tmp = input_pattern + "_tmp"
        shell(f"plink1.9 --keep-allele-order --allow-no-sex --bfile {input_pattern} "
              f"--update-chr {input.chrom} --make-bed --out {tmp}_1")
        if wildcards.dataset == "TwinGene":
            shell(f"plink1.9 --keep-allele-order --allow-no-sex --bfile {tmp}_1 --update-map {input.pos} --make-bed --out {output_pattern}")
        else:
            shell(f"plink1.9 --keep-allele-order --allow-no-sex --bfile {tmp}_1 --update-map {input.pos} --make-bed --out {tmp}_2")
            shell(f"plink1.9 --keep-allele-order --allow-no-sex --bfile {tmp}_2 --extract {input.pos} --make-bed --out {output_pattern}")
        shell(f"rm {tmp}*")


rule get_pos_chrom:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        chrom = "run_folder/update_files/{dataset}.chrom",
        pos   = "run_folder/update_files/{dataset}.pos"
    shell:
        """
        Rscript bin/format_input/update_reference.R {input.manifest} {output.chrom} 3
        """
