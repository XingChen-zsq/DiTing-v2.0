# visualization of carbon, nitrogen and sulfur cycle
import sketch  
rule diagrammatic_drawing:
    input:
        abundance_table=os.path.join(config["OUT_DIR"], 'pathways_relative_abundance.tab')
    output:
        carbon=os.path.join(config["OUT_DIR"], 'carbon_cycle_sketch.png'),
        nitrogen=os.path.join(config["OUT_DIR"], 'nitrogen_cycle_sketch.png'),
        sulfur=os.path.join(config["OUT_DIR"], 'sulfur_cycle_sketch.png'),
        DMSP=os.path.join(config["OUT_DIR"], 'DMSP_cycle_sketch.png')
    run:
        os.makedirs(config["OUT_DIR"], exist_ok=True)
        sketch.sketch(input.abundance_table, config["OUT_DIR"])

# visualization of heatmap
import heatmap  
rule heatmap_drawing:
    input:
        abundance_table=os.path.join(config["OUT_DIR"], 'pathways_relative_abundance.tab')
    output:
        carbon_cycle=os.path.join(config["OUT_DIR"], 'carbon_cycle_heatmap.pdf'),
        nitrogen_cycle=os.path.join(config["OUT_DIR"], 'nitrogen_cycle_heatmap.pdf'),
        sulfur_cycle=os.path.join(config["OUT_DIR"], 'sulfur_cycle_heatmap.pdf'),
        DMSP_cycle=os.path.join(config["OUT_DIR"], 'other_cycle_heatmap.pdf')
    run:
        os.makedirs(config["OUT_DIR"], exist_ok=True)
        heatmap.heatmap(input.abundance_table, config["OUT_DIR"])

