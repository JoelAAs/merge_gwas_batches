rule checkVCF_single:
    params:
        fasta_ref = config["fasta_ref"]
    input:
        "run_folder/single_vcf/{dataset}.vcf"
    output:
          "results/QC/single_vcf/{dataset}.check.af",
          "results/QC/single_vcf/{dataset}.check.dup",
          "results/QC/single_vcf/{dataset}.check.geno",
          "results/QC/single_vcf/{dataset}.check.log",
          "results/QC/single_vcf/{dataset}.check.mono",
          "results/QC/single_vcf/{dataset}.check.nonSnp",
          "results/QC/single_vcf/{dataset}.check.ref"
    log:
        "logs/QC/VCF_check_{dataset}.txt"
    shell:
        """
        python2.7 bin/vcf_handling/checkVCF.py -r {params.fasta_ref} -o results/QC/single_vcf/{wildcards.dataset} {input} &> {log}
        """

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
    log:
        "logs/QC/VCF_check_merged.txt"
    shell:
        """
        python2.7 bin/vcf_handling/checkVCF.py -r {params.fasta_ref} -o results/QC/merged {input}  &> {log}
        """