# kegg annotation
import csv
import logging
rule kegg_annotation:
    input:
        faa=os.path.join(config["CDHIT_DIR"], "{basename}_re.faa"),
        ko_list=os.path.join(config["KODB_DIR"], "ko_list"),
        profiles=os.path.join(config["KODB_DIR"], "profiles")
    output:
        os.path.join(config["KEGG_DIR"], "{basename}_annotation.tsv")
    params:
        cpu=config["cpu"]["kegg_annotaion"],
        temp_dir=lambda wildcards: f"./temp_{wildcards.basename}"  
    conda: "./metagenome.yaml"
    shell:
        """
        mkdir -p {params.temp_dir}
        exec_annotation -f detail-tsv -o {output} --cpu {params.cpu} -k {input.ko_list} -p {input.profiles} {input.faa} --tmp {params.temp_dir}
        rm -rf {params.temp_dir}
        """

# merge_and_process_annotationsg
rule merge_and_process_annotations:
    input:
        expand(os.path.join(config["KEGG_DIR"], "{basename}_annotation.tsv"), basename=config["BASENAMES"])
    output:
        os.path.join(config["KEGG_DIR"], "ko_merged.txt")
    run:
        process_and_merge_annotations(input, output[0])

def process_and_merge_annotations(input_files, output_file):
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(["#sample", "gene_id", "k_number"])

        for input_file in input_files:
            basename = os.path.basename(input_file).split('_annotation.tsv')[0]
            logging.info(f"handle {basename}:{input_file}")

            highest_score_rows = {}

            try:
                with open(input_file, 'r') as infile:
                    reader = csv.reader(infile, delimiter='\t')
                    headers = next(reader, None)
                    separator = next(reader, None)
                    if not headers or not separator:
                        logging.warning(f"{input_file} incorrectly formatted, skip...")
                        continue

                    for row in reader:
                        try:
                            gene_name, ko, threshold, score, e_value, ko_def = (
                                row[1], row[2], row[3], float(row[4]), float(row[5]), row[6]
                            )
                        except (IndexError, ValueError):
                            logging.warning(f"skip {input_file} incorrectly formatted: {row}")
                            continue

                        if e_value > 1e-05:
                            continue

                        if gene_name not in highest_score_rows or score > highest_score_rows[gene_name][2]:
                            highest_score_rows[gene_name] = [gene_name, ko, score]

            except Exception as e:
                logging.error(f"error reading file {input_file}: {e}")
                continue

            for gene_id, (gene_name, ko, score) in highest_score_rows.items():
                writer.writerow([basename, gene_name, ko])
                
    logging.info(f"The merge result has been written in {output_file}")
