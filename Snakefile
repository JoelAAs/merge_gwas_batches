configfile: "config.yml"
import re

include: "bin/format_input/formats.smk"
include: "bin/update_reference/update.smk"
include: "bin/merge/merge.smk"
include: "bin/vcf_handling/format_vcf.smk"
include: "bin/vcf_handling/check_vcf.smk"
include: "bin/qc_checks/Check_af.smk"

################ CONFIGURATION ################
wildcard_constraints:
    dataset = "[[a-zA-Z0-9-_\.]+"

################ RULES ################

rule all:
    input:
        "results/QC/merged.check.ref",
        "results/QC/AF_Gnomadcheck.csv"
        # expand("run_folder/single_vcf/{dataset}.vcf", dataset=config["datasets"]),
        # expand("results/QC/single_vcf/{dataset}.check.af", dataset=config["datasets"])
