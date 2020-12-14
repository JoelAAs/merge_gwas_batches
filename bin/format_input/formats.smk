## FLip files
rule get_flip_plus_to_forward:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        flipfile = "run_folder/flip/{dataset}_FORWARD.flip"
    shell:
        """
        Rscript bin/format_input/update_reference.R {input.manifest} {output.flipfile} 2
        """

rule get_flip_plus_from_top_bot:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        flipfile = "run_folder/flip/{dataset}_PLUS.flip"
    shell:
        """
        Rscript bin/format_input/update_reference.R {input.manifest} {output.flipfile} 0
        """


rule get_flip_top_bot_from_top:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        flip_top_bot = "run_folder/flip/{dataset}_TP.flip"
    shell:
        """
        Rscript bin/format_input/update_reference.R {input.manifest} {output.flip_top_bot} 1
        """


rule get_allele_order:
    input:
        manifest = lambda wc: config["manifest"][wc.dataset]
    output:
        allele_order = "run_folder/allele_order/{dataset}_{suffix}_order.txt"
    shell:
        """
        Rscript bin/format_input/get_allele_order.R {input.manifest} {wildcards.suffix} {output.allele_order} 
        """


## Format rules
rule set_bed_forward:
    """
    Update chromosome and position of ids given supplied dbsnp.map
    """
    input:
        bed = "run_folder/bed/{dataset}_PLUS.bed",
        flipfile = "run_folder/flip/{dataset}_FORWARD.flip"
    log:
        "logs/{dataset}_update_ref.log"
    output:
        bed = "run_folder/bed/{dataset}_FORWARD.bed",
        fam = "run_folder/bed/{dataset}_FORWARD.fam",
        bim = "run_folder/bed/{dataset}_FORWARD.bim"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --real-ref-alleles --allow-no-sex --bfile {input_pattern} --flip {input.flipfile} "
              f" --not-chr 0 --set-hh-missing --make-bed --out {output_pattern} &> {log}")

rule set_bed_plus:
    """
    Update chromosome and position of ids given supplied dbsnp.map
    """
    input:
        bed = "run_folder/bed/{dataset}_TP.bed",
        flipfile = "run_folder/flip/{dataset}_PLUS.flip"
    log:
        "logs/{dataset}_update_ref.log"
    output:
        bed = "run_folder/bed/{dataset}_PLUS.bed"
    run:
        input_pattern = re.sub("\\.bed", "", input[0])
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --real-ref-alleles --allow-no-sex --bfile {input_pattern} --flip {input.flipfile} "
              f"--make-bed --out {output_pattern} &> {log}")


rule set_bed_top_bot:
    input:
        bed          = "run_folder/bed/{dataset}_TOP.bed",
        flip_top_bot = "run_folder/flip/{dataset}_TP.flip"
    output:
        bed = "run_folder/bed/{dataset}_TP.bed"
    run:
        input_pattern = re.sub("\\.bed", "", input.bed)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --real-ref-alleles --allow-no-sex --bfile {input_pattern} "
              f"--flip {input.flip_top_bot} --make-bed --out {output_pattern}")


rule format_bed:
    """
    Format pad to bed
    """
    params:
        plink     = config["plink"]
    wildcard_constraints:
        suffix  = "(TOP|PLUS)",
    input:
        ped = "input/{dataset}_{suffix}.ped"
        # allele_order = "run_folder/allele_order/{dataset}_{suffix}_order.txt"
    output:
        bed = "run_folder/bed/{dataset}_{suffix}.bed"
    run:
        input_pattern = re.sub("\\.ped", "", input.ped)
        output_pattern = re.sub("\\.bed", "", output.bed)
        shell(f"plink1.9 --allow-no-sex --file {input_pattern} "
              f"--not-chr 0 --set-hh-missing "
              f"--make-bed --out {output_pattern}")
              # f"--a2-allele {input.allele_order} "
