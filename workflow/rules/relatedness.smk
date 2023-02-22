

rule vcftools_relatedness2:
    input:
        "results/merged_calls/all.vcf.gz",
    output:
        "results/relatedness/all.vcf.relatedness2"
    log:
        "results/logs/relatedness/vcftools_relatedness2.log"
    params:
        extra=config["vcftools_relatedness2"]["params"],
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools "
        "--gzvcf {input} "
        "--out results/relatedness/all.vcf "
        "{params.extra} "
        "--relatedness2 "
        "1>{log} 2>&1"


rule vcf_clone_detect:
    input:
        "results/merged_calls/all.vcf.gz",
    output:
        "results/relatedness/all.vcf.allelic_similarity.csv"
    log:
        "results/logs/relatedness/vcftools_relatedness2.log"
    params:
        threshold=config["vcf_clone_detect"]["threshold"],
        extra=config["vcf_clone_detect"]["params"],
    conda:
        "../envs/vcf_clone_detect.yaml"
    shell:
        "python workflow/scripts/vcf_clone_detect.py "
        "--vcf <(gunzip -c {input}) "
        "--output {output} "
        "--threshold {params.threshold} "
        "{params.extra} "
        "1>{log} 2>&1"


rule write_bam_list:
    input:
        expand("results/mapped/{sample}.bam", sample=samples.sample_id)
    output:
        files="results/relatedness/bam.filelist",
        labels="results/relatedness/bam.filelist.labels",
    log:
        "results/logs/relatedness/write_bam_list.log",
    shell:'''
        rm -fr {output.files} {output.labels}
        for f in {input};
        do
          echo $f >> {output.files}
          echo $f | sed -e 's@.*/@@' -e 's/.bam//' >> {output.labels}
        done
        '''


rule install_relatedness_programs:
    output:
        ".snakemake/install_relatedness_programs.done"
    log:
        "results/logs/relatedness/install_relatedness_programs.log",
    conda:
        "../envs/relatedness.yaml"
    threads: config["install_relatedness_programs"]["threads"]
    shell:'''
    (
    pwd=$PWD
    base=$(which python | sed -e 's@/bin/python@@')
    echo "base dir = $base"
    cd $base/lib
    rm -fr htslib ngsRelate pcangsd
    
    echo "##"
    echo "## Installing ngsRelate"
    echo "##"
    git clone --recursive https://github.com/SAMtools/htslib
    cd htslib/
    make -j {threads}
    cd ../
    git clone https://github.com/ANGSD/ngsRelate
    cd ngsRelate
    make -j {threads} HTSSRC=../htslib/
    cp ngsRelate ../../bin/
    cd ../
    
    echo "##"
    echo "## Installing PCAngsd"
    echo "##"
    git clone https://github.com/Rosemeis/pcangsd.git
    cd pcangsd
    python setup.py build_ext --inplace
    pip3 install -e .
    
    cd $pwd
    touch {output}
    ) 1>{log} 2>&1
    '''


rule angsd_for_NgsRelate:
    input:
        "results/relatedness/bam.filelist",
    output:
        mafs="results/relatedness/NgsRelate.angsd.mafs.gz",
        glf="results/relatedness/NgsRelate.angsd.glf.gz",
    log:
        "results/logs/relatedness/NgsRelate.angsd.log",
    params:
        extra=config["angsd_for_NgsRelate"]["params"],
    threads: config["angsd_for_NgsRelate"]["threads"]
    conda:
        "../envs/angsd.yaml"
    shell:
        "(angsd -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -nThreads {threads} {params.extra} -bam {input} -out results/relatedness/NgsRelate.angsd) 1>{log} 2>&1"


rule NgsRelate:
    input:
        mafs="results/relatedness/NgsRelate.angsd.mafs.gz",
        glf="results/relatedness/NgsRelate.angsd.glf.gz",
        labels="results/relatedness/bam.filelist.labels",
        programs=".snakemake/install_relatedness_programs.done"
    output:
        "results/relatedness/NgsRelate.results.tsv",
    log:
        "results/logs/relatedness/NgsRelate.log",
    params:
        extra=config["NgsRelate"]["params"],
    threads: config["NgsRelate"]["threads"]
    conda:
        "../envs/relatedness.yaml"
    shell:
        "(zcat {input.mafs} | cut -f5 |sed 1d > {input.mafs}.freq; "
        "N=$(cat {input.labels} | wc -l); "
        "ngsRelate -g {input.glf} -n $N -f {input.mafs}.freq -z {input.labels} -O {output} -p {threads} {params.extra}) 1>{log} 2>&1"


rule angsd_for_PCAngsd:
    input:
        "results/relatedness/bam.filelist",
    output:
        mafs="results/relatedness/PCAngsd.angsd.mafs.gz",
        beagle="results/relatedness/PCAngsd.angsd.beagle.gz",
    log:
        "results/logs/relatedness/PCAngsd.angsd.log",
    params:
        extra=config["angsd_for_PCAngsd"]["params"],
    threads: config["angsd_for_PCAngsd"]["threads"]
    conda:
        "../envs/angsd.yaml"
    shell:
        "(angsd -GL 2 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -nThreads {threads} {params.extra} -bam {input} -out results/relatedness/PCAngsd.angsd) 1>{log} 2>&1"


rule PCAngsd_IndAlleleFreq:
    input:
        beagle="results/relatedness/PCAngsd.angsd.beagle.gz",
        programs=".snakemake/install_relatedness_programs.done"
    output:
        "results/relatedness/PCAngsd.IndAlleleFreq.cov",
    log:
        "results/logs/relatedness/PCAngsd.IndAlleleFreq.log",
    params:
        extra=config["PCAngsd_IndAlleleFreq"]["params"],
    threads: config["PCAngsd_IndAlleleFreq"]["threads"]
    conda:
        "../envs/relatedness.yaml"
    shell:
        "(pcangsd --threads {threads} {params.extra} --beagle {input.beagle} --out results/relatedness/PCAngsd.IndAlleleFreq) 1>{log} 2>&1"


rule PCAngsd_WithOutIndAlleleFreq:
    input:
        beagle="results/relatedness/PCAngsd.angsd.beagle.gz",
        programs=".snakemake/install_relatedness_programs.done"
    output:
        "results/relatedness/PCAngsd.WithOutIndAlleleFreq.cov",
    log:
        "results/logs/relatedness/PCAngsd.WithOutIndAlleleFreq.log",
    params:
        extra=config["PCAngsd_WithOutIndAlleleFreq"]["params"],
    threads: config["PCAngsd_WithOutIndAlleleFreq"]["threads"]
    conda:
        "../envs/relatedness.yaml"
    shell:
        "(pcangsd --threads {threads} {params.extra} --beagle {input.beagle} --out results/relatedness/PCAngsd.WithOutIndAlleleFreq --iter 0) 1>{log} 2>&1"


rule PCAngsd_Admixture:
    input:
        beagle="results/relatedness/PCAngsd.angsd.beagle.gz",
        programs=".snakemake/install_relatedness_programs.done"
    output:
        "results/relatedness/PCAngsd.Admixture.cov",
    log:
        "results/logs/relatedness/PCAngsd.Admixture.log",
    params:
        extra=config["PCAngsd_Admixture"]["params"],
    threads: config["PCAngsd_Admixture"]["threads"]
    conda:
        "../envs/relatedness.yaml"
    shell:
        "(pcangsd --threads {threads} {params.extra} --beagle {input.beagle} --out results/relatedness/PCAngsd.Admixture --admix --admix_alpha 50) 1>{log} 2>&1"

