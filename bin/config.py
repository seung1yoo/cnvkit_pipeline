
import os
import sys

class Utils:
    def __init__(self, config):
        self.config = config

    def get_targets(self, targets):
        ls = list()
        if 'symlink_bam' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/inbam/{sample}/{sample}.bam')
                ls.append(f'analysis/inbam/{sample}/{sample}.bai')
        if 'run_qualimap_bamqc' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/qualimap/{sample}/qualimapReport.html')
        if 'run_qualimap_multibamqc' in targets:
            ls.append(f'analysis/qualimap/multi/input_data.txt')
            ls.append(f'analysis/qualimap/multi/multisampleBamQcReport.html')
        if 'bedtools_coverage' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/bedtools/{sample}/coverage.txt')
        if 'bedtools_coverage_bed' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/bedtools/{sample}/coverage.bed.txt')
        if 'bedtools_genomecov' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/bedtools/{sample}/genomecov.txt')
        if 'samtools_depth_genome' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/samtools/{sample}/depth.genome.txt')
        if 'samtools_depth_target' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/samtools/{sample}/depth.target.txt')
        if 'covtobed' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/covtobed/{sample}/covered.10x.bed')
        if 'covtobed_allsamples' in targets:
            ls.append(f'analysis/covtobed/allsamples/covered.10x.bed')
        if 'bedtools_merge' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/covtobed/{sample}/covered.10x.merge.bed')
        if 'bedtools_intersect' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/covtobed/{sample}/covered.10x.target.bed')
        if 'bedtools_intersect_co' in targets:
            ls.append(f'analysis/covtobed/{sample}/covered.10x.target.bed')
        if 'cnvkit_batch' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}.cnn')
                ls.append(f'analysis/cnvkit/{sample}/{sample}.cnr')
                ls.append(f'analysis/cnvkit/{sample}/{sample}.cns')
        if 'cnvkit_segment' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}.cns')
        if 'cnvkit_call' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}.call.cns')
        if 'cnvkit_call_bed' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}.call.bed')
        if 'cnvkit_call_anno' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}.call.anno.bed')
        if 'cnvkit_gainloss' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}.gene.gainloss.txt')
        if 'cnvkit_scatter' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}-scatter.pdf')
                ls.append(f'analysis/cnvkit/{sample}/{sample}-scatter.png')
        if 'cnvkit_diagram' in targets:
            for sample in self.config['ordered_samples']:
                ls.append(f'analysis/cnvkit/{sample}/{sample}-diagram.pdf')
                ls.append(f'analysis/cnvkit/{sample}/{sample}-diagram.png')
        if 'parse_cnvkit' in targets:
            ls.append(f'analysis/cnvkit/compbed/ref.ipsc57LBvsRef.bed')
            ls.append(f'analysis/cnvkit/compbed/comp.ipsc57LBvsRef.bed')
        if 'cnvkit_heatmap' in targets:
            ls.append(f'analysis/cnvkit/heatmap.png')
        return ls



class MakeConfig:
    def __init__(self, id_dic, dirPath_bam, workdir):
        self.id_dic = id_dic
        self.init_config_dic()
        self.addpath_bam(dirPath_bam)
        self.workdir = workdir
        #
        self.write_config('config.yaml')


    def init_config_dic(self):
        self.config_dic = dict()
        for tbi_id, cst_id in self.id_dic.items():
            self.config_dic.setdefault('samples', {}).setdefault(cst_id, {}).setdefault('tbi_id', tbi_id)
            self.config_dic.setdefault('samples', {}).setdefault(cst_id, {}).setdefault('bam', None)
            self.config_dic.setdefault('samples', {}).setdefault(cst_id, {}).setdefault('bai', None)

    def addpath_bam(self, _path):
        for tbi_id, cst_id in self.id_dic.items():
            #bam = os.path.join(_path, tbi_id, f'{tbi_id}.dedup.bam')
            #bai = os.path.join(_path, tbi_id, f'{tbi_id}.dedup.bai')
            bam = os.path.join(_path, tbi_id, f'{tbi_id}.merge.bam')
            bai = os.path.join(_path, tbi_id, f'{tbi_id}.merge.bam.bai')
            if os.path.isfile(bam) and os.path.isfile(bai):
                self.config_dic['samples'][cst_id]['bam'] = os.path.abspath(bam)
                self.config_dic['samples'][cst_id]['bai'] = os.path.abspath(bai)
            else:
                print('Could not found : ', bam)
                print('or')
                print('Could not found : ', bai)
                sys.exit()
        #

    def write_config(self, outfn):
        outfh = open(outfn, 'w')
        outfh.write('\n')
        outfh.write('workdir: {0}\n'.format(self.workdir))
        outfh.write('target_bed: /BiO/BioResources/References/Human/hg19/targetkit/SureSelect_Human_All_Exon_V5.bed\n')
        outfh.write('genome_fasta: /BiO/BioResources/References/Human/hg19/hg19.fa\n')
        outfh.write('access_bed: /BiO/BioTools/cnvkit/data/access-5k-mappable.hg19.bed\n')
        outfh.write('ref_flat: /BiO/BioPeople/baekip/BioResource/hg19/refFlat.txt\n')
        outfh.write('ipsc57_bed: bin/ipsc_57_LBvsRef_Region.bed\n')
        outfh.write('\n')
        outfh.write('ordered_samples:\n')
        for sample_id, info_dic in sorted(self.config_dic['samples'].items()):
            outfh.write('  - {0}\n'.format(sample_id))
        outfh.write('\n')
        outfh.write('samples:\n')
        for sample_id, info_dic in self.config_dic['samples'].items():
            outfh.write('  {0}:\n'.format(sample_id))
            outfh.write('    tbi_id: {0}\n'.format(info_dic['tbi_id']))
            outfh.write('    bam: {0}\n'.format(info_dic['bam']))
            outfh.write('    bai: {0}\n'.format(info_dic['bai']))
        outfh.close()




def main():
    id_dic = {"TN1808L0024-10":"H9ESP_44-10",
              "TN1808L0024-1":"H9ESP_44-1",
              "TN1808L0024-2":"H9ESP_44-2",
              "TN1808L0024-3":"H9ESP_44-3",
              "TN1808L0024-4":"H9ESP_44-4",
              "TN1808L0024-5":"H9ESP_44-5",
              "TN1808L0024-6":"H9ESP_44-6",
              "TN1808L0024-7":"H9ESP_44-7",
              "TN1808L0024-8":"H9ESP_44-8",
              "TN1808L0024-9":"H9ESP_44-9",
              "TN1808L0025-10":"H9ESP_75-10",
              "TN1808L0025-1":"H9ESP_75-1",
              "TN1808L0025-2":"H9ESP_75-2",
              "TN1808L0025-3":"H9ESP_75-3",
              "TN1808L0025-4":"H9ESP_75-4",
              "TN1808L0025-5":"H9ESP_75-5",
              "TN1808L0025-6":"H9ESP_75-6",
              "TN1808L0025-7":"H9ESP_75-7",
              "TN1808L0025-8":"H9ESP_75-8",
              "TN1808L0025-9":"H9ESP_75-9",
              "TN1808L0030-10":"HPS0076_48-10",
              "TN1808L0030-1":"HPS0076_48-1",
              "TN1808L0030-2":"HPS0076_48-2",
              "TN1808L0030-3":"HPS0076_48-3",
              "TN1808L0030-4":"HPS0076_48-4",
              "TN1808L0030-5":"HPS0076_48-5",
              "TN1808L0030-6":"HPS0076_48-6",
              "TN1808L0030-7":"HPS0076_48-7",
              "TN1808L0030-8":"HPS0076_48-8",
              "TN1808L0030-9":"HPS0076_48-9",
              "TN1808L0031-10":"HPS0076_68-10",
              "TN1808L0031-1":"HPS0076_68-1",
              "TN1808L0031-2":"HPS0076_68-2",
              "TN1808L0031-3":"HPS0076_68-3",
              "TN1808L0031-4":"HPS0076_68-4",
              "TN1808L0031-5":"HPS0076_68-5",
              "TN1808L0031-6":"HPS0076_68-6",
              "TN1808L0031-7":"HPS0076_68-7",
              "TN1808L0031-8":"HPS0076_68-8",
              "TN1808L0031-9":"HPS0076_68-9"}
    #dirPath_bam = "/BiO/BioPeople/baekip/BioProjects/StemCell_Project/MFDS-Human-C1-2018-08-TBD171103/result/07-1_picard_dedup"
    dirPath_bam = "/BiO/BioPeople/baekip/BioProjects/StemCell_Project/MFDS-Human-C1-2018-08-TBD171103/result/06-1_picard_merge"
    workdir = "/BiO/BioPeople/baekip/BioProjects/StemCell_Project/MFDS-Human-C1-2018-08-TBD171103/202012301118_cnvkit"
    #
    config_obj = MakeConfig(id_dic, dirPath_bam, workdir)




if __name__=='__main__':
    main()


