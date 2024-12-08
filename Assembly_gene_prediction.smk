# MEGAHIT
rule megahit:
    input:
        read1=os.path.join(config["READS_DIR"], "{basename}_1" + config["READS_SUF"]),
        read2=os.path.join(config["READS_DIR"], "{basename}_2" + config["READS_SUF"])
    output:
        contigs=os.path.join(config["ASSEMBLY_DIR"], "{basename}.fa")
    conda:
        "./metagenome.yaml"
    log:
        os.path.join(config["ASSEMBLY_DIR"], "{basename}.log")  
    threads: config["threads"]["megahit"]
    shell:
        """
        megahit -1 {input.read1} -2 {input.read2} -o {output.contigs}.tmp -t {threads} &> {log} \
        && mv {output.contigs}.tmp/final.contigs.fa {output.contigs} \
        && rm -r {output.contigs}.tmp
        """
# Prodigal
rule prodigal:
    input:
        contigs=os.path.join(config["ASSEMBLY_DIR"], "{basename}.fa")
    output:
        proteins=os.path.join(config["PRODIGAL_DIR"], "{basename}.faa"),  
        genes=os.path.join(config["PRODIGAL_DIR"], "{basename}.ffn"),  
        gbk=os.path.join(config["PRODIGAL_DIR"], "{basename}.gbk")
    params:
        mode="meta"  
    conda: "./metagenome.yaml"
    log:
        os.path.join(config["PRODIGAL_DIR"], "{basename}.log")
    threads: config["threads"]["prodigal"]
    shell:
        """
        prodigal -i {input.contigs} -a {output.proteins} -d {output.genes} -o {output.gbk} -p {params.mode} -t {threads} &> {log}
        """

# CD-HIT 
rule cdhit:
    input:
        proteins=os.path.join(config["PRODIGAL_DIR"], "{basename}.faa")  
    output:
        non_redundant=os.path.join(config["CDHIT_DIR"], "{basename}_re.faa")
    params:
        identity_threshold=0.95,
        threads=config["threads"]["cdhit"]   
    conda: "./metagenome.yaml"
    log:
        os.path.join(config["CDHIT_DIR"], "{basename}.log")
    shell:
        """
        cd-hit -i {input.proteins} -o {output.non_redundant} -c {params.identity_threshold} -n 5 -M 40000 -d 0 -T {params.threads} &> {log}
        """

# Filter_nucleotides
rule filter_nucleotides:
    input:
        re_faa=os.path.join(config["CDHIT_DIR"], "{basename}_re.faa"),  
        ffn=os.path.join(config["PRODIGAL_DIR"], "{basename}.ffn")  
    output:
        re_ffn=os.path.join(config["FILTERED_DIR"], "{basename}_filtered.ffn")  
    run:
        filter_nucleotides_by_protein(input.re_faa, input.ffn, output.re_ffn)

def filter_nucleotides_by_protein(re_faa, ffn, re_ffn):
    kept_sequences = set()
    with open(re_faa, 'r') as faa_file:
        for line in faa_file:
            if line.startswith('>'):
                seq_name = line.split()[0].lstrip('>')
                kept_sequences.add(seq_name)
    with open(ffn, 'r') as ffn_file, open(re_ffn, 'w') as out_file:
        write_flag = False
        for line in ffn_file:
            if line.startswith('>'):
                seq_name = line.split()[0].lstrip('>')
                write_flag = seq_name in kept_sequences
            if write_flag:
                out_file.write(line)
