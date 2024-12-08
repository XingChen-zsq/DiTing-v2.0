# merge gene relative abundance table with gene kegg annotation table
import logging
rule merge_abun_ko:
    input:
        abun_files=expand(os.path.join(config["GENE_ABUN_DIR"], "{basename}.abundance"), basename=config["BASENAMES"]),
        ko_merged=os.path.join(config["KEGG_DIR"], 'ko_merged.txt')
    output:
        ko_abun=os.path.join(config["KEGG_DIR"], 'ko_abun.txt')
    run:
        logging.info("Merge KEGG annotations and gene relative abundances")
        merge_abun_ko(config["GENE_ABUN_DIR"], input.ko_merged, output.ko_abun)

def merge_abun_ko(GENE_ABUN_DIR, ko_merged, output):
    logging.info('Merge abundance table with KEGG table')
    with open(output, 'w') as fo:
        fo.write('#sample\tk_number\trelative_abundance\tgene_id\n')
    
    abun_tab_dict = {}
    for abun in os.listdir(GENE_ABUN_DIR):
        if abun.endswith('.abundance'):
            basename = abun.rsplit('.', 1)[0]  
            abun_path = os.path.join(GENE_ABUN_DIR, abun)
            with open(abun_path) as fi:
                next(fi)  
                for line in fi:
                    line = line.strip('\n')
                    gene_id, abundance = line.split('\t')
                    key = str(basename) + '+' + str(gene_id)  
                    abun_tab_dict[key] = abundance  
    a = []
    with open(ko_merged, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                a.append(line.strip())
    for item in sorted(a):
        basename = item.split('\t')[0]
        gene_id = item.split('\t')[1]
        k_number = item.split('\t')[2]
        key = basename + '+' + gene_id
        abundance = abun_tab_dict.get(key, "NA")  
        with open(output, 'a') as fo:
            fo.write(basename + '\t' + k_number + '\t' + abundance + '\t' + gene_id + '\n')


# Produce table of ko abundance among samples
rule table_of_ko_abundance_among_samples:
    input:
        ko_abun_txt=os.path.join(config["KEGG_DIR"], 'ko_abun.txt')
    output:
        ko_abun_among_samples=os.path.join(config["KEGG_DIR"], 'ko_abundance_among_samples.tab')
    run:
        table_of_ko_abundance_among_samples(input.ko_abun_txt, output.ko_abun_among_samples)

def table_of_ko_abundance_among_samples(ko_abun_txt, output):
    sampleKnumber_to_abundance = {}
    samples = []
    k_numbers = []
    logging.info('Produce table of ko abundance among samples')
    
    with open(ko_abun_txt) as fi:
        for line in fi:
            line = line.strip()   
            if line.startswith('#'):
                continue
            else:
                sample = line.split('\t')[0]
                k_number = line.split('\t')[1]
                abundance = line.split('\t')[2]
                key = sample + '+' + k_number
                if key in sampleKnumber_to_abundance:
                    sampleKnumber_to_abundance[key] += float(abundance)
                else:
                    sampleKnumber_to_abundance[key] = float(abundance)
                samples.append(sample)
                k_numbers.append(k_number)

    samples2 = list(set(samples))
    samples2.sort()
    k_numbers2 = list(set(k_numbers))   
    k_numbers2.sort()
   
    with open(output, 'w') as fo:
        fo.write('k_number')  
        for sample in samples2:
            fo.write('\t' + sample)  
        fo.write('\n')

    for k_number in k_numbers2:
        with open(output, 'a') as fo:
            fo.write(k_number)  
            for sample in samples2:
                key = sample + '+' + k_number
                if key not in sampleKnumber_to_abundance:
                    sampleKnumber_to_abundance[key] = float(0)  
                fo.write('\t' + str(sampleKnumber_to_abundance[key]))  
            fo.write('\n')


# Generate hierarchical table abundance among samples
rule hierarchical_ko_abundance:
    input:
        ko_abundance_among_samples=os.path.join(config["KEGG_DIR"], 'ko_abundance_among_samples.tab'),
        KO_affilated_to_biogeochemical_cycle=os.path.join(config["TABLE"], 'KO_affilated_to_biogeochemical_cycle.tab')
    output:
        pathways_relative_abundance_gene_level=os.path.join(config["KEGG_DIR"], 'pathways_relative_abundance_gene_level.tab')
    run:
        logging.info("Generate hierarchical table")
        hierarchical_ko_abundance_among_samples(input.ko_abundance_among_samples,
                                                 input.KO_affilated_to_biogeochemical_cycle,
                                                 output.pathways_relative_abundance_gene_level)

def hierarchical_ko_abundance_among_samples(ko_abundance_among_samples, KO_affilated_to_biogeochemical_cycle, output):
    Knumber_to_abundance = {}
    with open(ko_abundance_among_samples) as fi:
        for line in fi:
            line = line.rstrip()
            lines = line.split('\t')
            k_number = lines[0]
            lines.pop(0)
            value = '\t'.join(lines)
            Knumber_to_abundance[k_number] = value

    with open(output, 'w') as fo:
        fo.write('') 

    with open(KO_affilated_to_biogeochemical_cycle) as fi:
        for line in fi:
            line = line.rstrip()
            k_number = line.split('\t')[2]
            if k_number not in Knumber_to_abundance:
                with open(output, 'a') as fo:
                    fo.write(line + '\n')
            else:
                with open(output, 'a') as fo:
                    fo.write(line + '\t' + str(Knumber_to_abundance[k_number]) + '\n')


# Calculate relative abundance of pathways
import pathway  
rule calculate_pathway_abundance:
    input:
        ko_abun=os.path.join(config["KEGG_DIR"], 'ko_abun.txt')
    output:
        pathway_abundance=os.path.join(config["OUT_DIR"], 'pathways_relative_abundance.tab')
    run:
        logging.info("Calculate abundance of pathways")
        function_order = pathway.function_order  
        pathway_parser(input.ko_abun, output.pathway_abundance, function_order)

def pathway_parser(ko_abun, pathway_abundance, function_order):
    relative_abundance = {}
    genome_data = []
    with open(ko_abun, "r") as f:
        for line in f:
            line = line.rstrip()
            info = line.split()
            if line.startswith('#'):
                continue
            else:
                line_key = info[0] + '_' + info[1]
                relative_abundance[line_key] = relative_abundance.get(line_key, 0) + float(info[2])
                genome_data.append(info[0])
    genome_data = list(dict.fromkeys(genome_data))
    data_to_transpose = []
    header = ['#Pathway'] + function_order  
    data_to_transpose.append(header)
    for k in genome_data:
        pathway_abundant = pathway.Pathway(relative_abundance, k)  
        pathway_data = pathway_abundant.solve_pathway()
        out_list = [k] + [pathway_data.get(i, 0) for i in function_order]
        data_to_transpose.append(out_list)

    transposed_data = list(zip(*data_to_transpose))

    with open(pathway_abundance, "w") as out_file:
        for row in transposed_data:
            out_file.write("\t".join(map(str, row)) + "\n")
