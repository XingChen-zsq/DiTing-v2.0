import os
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import shutil

def heatmap(abundance_table, OUT_DIR):
    os.makedirs('heatmap_tmp', exist_ok=True)
    dict_table = {}
    head = ''
    with open(abundance_table) as table:
        for line in table:
            line = line.strip('\n')
            if line.startswith('#'):
                head = line
            else:
                pathway = line.split('\t')[0]
                abundance = line.split('\t')
                del abundance[0]
                values = '\t'.join(abundance)
                dict_table[pathway] = values

    # 定义路径循环
    cycles = {
        'carbon_cycle': ['Photosystem II (psbABCDEF)', 'Photosystem I (psaABCDEF)', 'Cytochrome b6/f complex (petABCDGLMN)', 
                         'Anoxygenic photosystem II (pufML)', 'Anoxygenic photosystem I (pscABCD)', 'RuBisCo'],
        'nitrogen_cycle': ['Dissimilatory nitrate reduction, nitrate -> nitrite (narGHI or napAB)'],
        'sulfur_cycle': ['Assimilatory sulfate reduction, sulfate -> sulfite', 'Assimilatory sulfate reduction, sulfite -> sulfide (cysJI or sir)'],
        'other_cycle': ['F-type ATPase', 'V/A-type ATPase']
    }

    # 根据路径写入数据
    for cycle_name, pathways in cycles.items():
        with open(f'heatmap_tmp/{cycle_name}.tab', 'w') as cycle_out:
            cycle_out.write(head + '\n')
            for pathway in pathways:
                cycle_out.write(pathway + '\t' + dict_table[pathway] + '\n')

    # 生成并保存热图
    for cycle_name in cycles.keys():
        plt.cla()
        file_in = f'heatmap_tmp/{cycle_name}.tab'
        data = pd.read_table(file_in, index_col=0)
        ax = sns.heatmap(data, cmap='coolwarm', xticklabels=True, yticklabels=True, square=True)
        fig = ax.get_figure()
        fig.set_size_inches(8.27, 11.69)
        fig.savefig(f"{OUT_DIR}/{cycle_name}_heatmap.pdf", bbox_inches='tight', dpi=600)
        plt.close()

    # 清理临时文件夹
    shutil.rmtree('heatmap_tmp')

