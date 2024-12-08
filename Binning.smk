from glob import glob
import logging
# metawrap_binning
rule metawrap_binning:
    input:
        read_1=os.path.join(config["READS_DIR"], "{basename}_1" + config["READS_SUF"]),
        read_2=os.path.join(config["READS_DIR"], "{basename}_2" + config["READS_SUF"]),
        contigs=os.path.join(config["ASSEMBLY_DIR"], "{basename}.fa")
    output:
        metabat_bins=directory(os.path.join(config["METAWRAP_DIR"], "{basename}_initial_bin/metabat2_bins/")),
        maxbin_bins=directory(os.path.join(config["METAWRAP_DIR"], "{basename}_initial_bin/maxbin2_bins/")),
        concoct_bins=directory(os.path.join(config["METAWRAP_DIR"], "{basename}_initial_bin/concoct_bins/"))
    params:
        mem=80  
    conda:
        "./genomic.yaml" 
    threads: config["threads"]["metawrap_binning"] 
    shell:
        """
        metawrap binning -o {config[METAWRAP_DIR]}/{wildcards.basename}_initial_bin -t {threads} -m {params.mem} -a {input.contigs} --metabat2 --maxbin2 --concoct {input.read_1} {input.read_2}
        """
# metawrap_bin_refinement
rule metawrap_bin_refinement:
    input:
        metabat_bins=os.path.join(config["METAWRAP_DIR"], "{basename}_initial_bin/metabat2_bins/"),
        maxbin_bins=os.path.join(config["METAWRAP_DIR"], "{basename}_initial_bin/maxbin2_bins/"),
        concoct_bins=os.path.join(config["METAWRAP_DIR"], "{basename}_initial_bin/concoct_bins/")
    output:
        refined_bins_dir=directory(os.path.join(config["METAWRAP_DIR"], "{basename}_REFINEMENT"))
    params:
        completeness=50,
        contamination=10,
        threads=10
    conda:
        "./genomic.yaml"
    threads: config["threads"]["metawrap_refinement"]
    shell:
        """
        metawrap bin_refinement -o {output.refined_bins_dir} \
                                -t {threads} \
                                -A {input.metabat_bins} \
                                -B {input.maxbin_bins} \
                                -C {input.concoct_bins} \
                                -c {params.completeness} \
                                -x {params.contamination}
        """

# Rename and copy the bin file 
rule rename_and_copy_bins:
    input:
        refinement_dir=os.path.join(config["METAWRAP_DIR"], "{basename}_REFINEMENT")
    output:
        directory(os.path.join(config["FINAL_BINS_DIR"], "{basename}"))  
    params:
        script="rename_and_copy_bins.py"
    shell:
        """
        python {params.script} {input.refinement_dir} {output} {wildcards.basename}
        """

rule gather_all_bins:
    input:
        expand(os.path.join(config["FINAL_BINS_DIR"], "{basename}"), basename=config["BASENAMES"])
    output:
        touch(os.path.join(config["ALL_BINS_DIR"], "completed.txt"))
    params:
        output_dir=config["ALL_BINS_DIR"]
    run:
        import shutil
        os.makedirs(params.output_dir, exist_ok=True)
        for bin_dir in input:
            for fa_file in glob(os.path.join(bin_dir, "*.fa")):
                dest_path = os.path.join(params.output_dir, os.path.basename(fa_file))
                shutil.copy(fa_file, dest_path)
                logging.info(f"Copied {fa_file} to {dest_path}")

# metawrap_quantify_bins
rule quantify_bins:
    input:
        bins_dir=config["ALL_BINS_DIR"],
        completed_txt=os.path.join(config["ALL_BINS_DIR"], "completed.txt"), 
        read1=glob(os.path.join(config["READS_DIR"], "*_1.fastq")),
        read2=glob(os.path.join(config["READS_DIR"], "*_2.fastq")),
        assemblies=glob(os.path.join(config["ASSEMBLY_DIR"], "*.fa"))
    output:
        abundance_file=os.path.join(config["BIN_ABUNDANCE_DIR"], "bin_abundance_table.tab")
    conda:
        "./genomic.yaml"
    threads: config["threads"]["quantify_bins"]
    shell:
        """
        metawrap quant_bins -b {input.bins_dir} -t {threads} -o {config[BIN_ABUNDANCE_DIR]} -a {input.assemblies} {input.read1} {input.read2}
        """
