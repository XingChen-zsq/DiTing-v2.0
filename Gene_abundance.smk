# BWA Index
rule bwa_index:
    input:
        fasta=os.path.join(config["FILTERED_DIR"], "{basename}_filtered.ffn")
    output:
        bwt=os.path.join(config["BWA_INDEX"], "{basename}.bwt"),
        pac=os.path.join(config["BWA_INDEX"], "{basename}.pac"),
        ann=os.path.join(config["BWA_INDEX"], "{basename}.ann"),
        amb=os.path.join(config["BWA_INDEX"], "{basename}.amb"),
        sa=os.path.join(config["BWA_INDEX"], "{basename}.sa")
    log:
        os.path.join(config["BWA_INDEX"], "{basename}.log")
    conda: "./metagenome.yaml"
    shell:
        """
        bwa index -p {config[BWA_INDEX]}/{wildcards.basename} {input.fasta} &> {log}
        """

# BWA Mapping
rule bwa_mem:
    input:
        index_files=os.path.join(config["BWA_INDEX"], "{basename}.bwt"),
        read1=os.path.join(config["READS_DIR"], "{basename}_1" + config["READS_SUF"]),
        read2=os.path.join(config["READS_DIR"], "{basename}_2" + config["READS_SUF"])
    output:
        sam=os.path.join(config["MAPPING"], "{basename}.sam")
    conda: "./metagenome.yaml"
    log:
        os.path.join(config["MAPPING"], "{basename}.bwa_mem.log")
    threads: config["threads"]["bwa"]
    shell:
        """
        bwa mem -t {threads} -v 2 {config[BWA_INDEX]}/{wildcards.basename} {input.reads} > {output.sam} 2> {log}
        """


# Pileup 
rule pileup:
    input:
        sam=os.path.join(config["MAPPING"], "{basename}.sam")
    output:
        output_file=os.path.join(config["PILEUP"], "{basename}.pileup") 
    conda: "./metagenome.yaml"
    log:
        os.path.join(config["PILEUP"], "{basename}.pileup.log")
    shell:
        """
        pileup.sh in={input.sam} out={output.output_file} &> {log}
        """

# Calculate gene abundance
rule gene_abundance:
    input:
        pileup=os.path.join(config["PILEUP"], "{basename}.pileup")
    output:
        abundance=os.path.join(config["GENE_ABUN_DIR"], "{basename}.abundance")
    run:
        gene_relative_abun(input.pileup, wildcards.basename, config["GENE_ABUN_DIR"])

def gene_relative_abun(pileup_file, basename, GENE_ABUN_DIR):
    total_ave_fold = 0.0
    file_out = os.path.join(GENE_ABUN_DIR, basename + '.abundance')
    with open(file_out, 'w') as fo:
        fo.write("#ID\tgene_abundance\n")
    with open(pileup_file, 'r') as fi:
        for line in fi:
            if line.startswith('#ID'):
                continue
            else:
                ave_fold = line.split('\t')[1]
                total_ave_fold += float(ave_fold)
    with open(pileup_file, 'r') as fi:
        for line in fi:
            if line.startswith('#ID'):
                continue
            else:
                gene_id = line.split('\t')[0]
                ave_fold = line.split('\t')[1]
                gene_abund = (float(ave_fold) / total_ave_fold) * 1000000
                with open(file_out, 'a') as fo:
                    fo.write(f"{gene_id}\t{gene_abund}\n")

