

rule append_maf_gnomad:
    params:
        gnomad = "/media/joel/Encrypted/gnomad.genomes.r2.1.1.sites.vcf.bgz"
    input:
          "results/QC/merged.check.af"
    output:
          "results/QC/AF_Gnomadcheck.csv"
    shell:
         """
         python bin/qc_checks/ExploreGnomad.py {params.gnomad} {input} {output}
         """
