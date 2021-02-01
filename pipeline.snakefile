
import os

from snakemake.utils import min_version
min_version("5.31.1")

from bin.config import Utils

configfile: "config.yaml"
utils = Utils(config)

targets = list()
targets.append('symlink_bam')
targets.append('run_qualimap_bamqc')
targets.append('run_qualimap_multibamqc')
#targets.append('bedtools_coverage')
#targets.append('bedtools_coverage_bed')
#targets.append('bedtools_genomecov')
#targets.append('samtools_depth_genome')
#targets.append('samtools_depth_target')
targets.append('covtobed')
#targets.append('covtobed_allsamples')
targets.append('bedtools_merge')
targets.append('bedtools_intersect')
targets.append('bedtools_intersect_co')
targets.append('cnvkit_batch')
#targets.append('cnvkit_segment')
targets.append('cnvkit_call')
targets.append('cnvkit_call_bed')
targets.append('cnvkit_call_anno')
targets.append('cnvkit_gainloss')
targets.append('cnvkit_scatter')
targets.append('cnvkit_diagram')
targets.append('cnvkit_heatmap')
#targets.append('parse_cnvkit')

rule all:
    input: utils.get_targets(targets)

rule symlink_bam:
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample]["bam"],
        bai = lambda wildcards: config["samples"][wildcards.sample]["bai"]
    output:
        bam = 'analysis/inbam/{sample}/{sample}.bam',
        bai = 'analysis/inbam/{sample}/{sample}.bai'
    shell:
        "ln -s {input.bam} {output.bam} && "
        "ln -s {input.bai} {output.bai}"

rule bedtools_coverage:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam'
    output:
        coverage = 'analysis/bedtools/{sample}/coverage.txt'
    params:
        bed = config['target_bed']
    shell:
        "bedtools coverage"
        " -abam {input.bam}"
        " -b {params.bed}"
        " > {output.coverage}"

rule bedtools_coverage_bed:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam'
    output:
        coverage = 'analysis/bedtools/{sample}/coverage.bed.txt'
    params:
        bed = config['target_bed']
    shell:
        "bedtools coverage"
        " -a {params.bed}"
        " -b {input.bam}"
        " > {output.coverage}"

rule bedtools_genomecov:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam'
    output:
        genomecov = 'analysis/bedtools/{sample}/genomecov.txt'
    shell:
        "bedtools genomecov"
        " -ibam {input.bam}"
        " -bg"
        " > {output.genomecov}"

rule samtools_depth_genome:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam'
    output:
        depth = 'analysis/samtools/{sample}/depth.genome.txt'
    shell:
        "samtools depth"
        " {input.bam}"
        " -o {output.depth}"

rule samtools_depth_target:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam'
    output:
        depth = 'analysis/samtools/{sample}/depth.target.txt'
    params:
        target = config['target_bed']
    shell:
        "samtools depth"
        " {input.bam}"
        " -b {params.target}"
        " -o {output.depth}"

rule covtobed:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam'
    output:
        bed = 'analysis/covtobed/{sample}/covered.10x.bed'
    params:
        min_cov = 10,
        max_cov = 100000,
        min_len = 1
    shell:
        "covtobed "
        " -m {params.min_cov}"
        " -x {params.max_cov}"
        " -l {params.min_len}"
        " {input.bam}"
        " > {output.bed}"

#rule covtobed_allsamples:
#    input:
#        bams = expand('analysis/inbam/{sample}/{sample}.bam', sample=config['ordered_samples'])
#    output:
#        bed = 'analysis/covtobed/allsamples/covered.10x.bed'
#    params:
#        min_cov = 10,
#        max_cov = 100000,
#        min_len = 1
#    shell:
#        "covtobed "
#        " -m {params.min_cov}"
#        " -x {params.max_cov}"
#        " -l {params.min_len}"
#        " {input.bams}"
#        " > {output.bed}"

rule bedtools_merge:
    input:
        bed = 'analysis/covtobed/{sample}/covered.10x.bed'
    output:
        bed = 'analysis/covtobed/{sample}/covered.10x.merge.bed'
    shell:
        "bedtools merge"
        " -i {input.bed}"
        " > {output.bed}"

rule bedtools_intersect:
    input:
        bed = 'analysis/covtobed/{sample}/covered.10x.merge.bed'
    output:
        bed = 'analysis/covtobed/{sample}/covered.10x.target.bed'
    params:
        target = config['target_bed']
    shell:
        "bedtools intersect"
        " -a {input.bed}"
        " -b {params.target}"
        " > {output.bed}"

rule bedtools_intersect_co:
    input:
        #beds = expand('analysis/covtobed/{sample}/covered.10x.merge.bed', sample=config['ordered_samples'])
        beds = expand('analysis/covtobed/{sample}/covered.10x.target.bed', sample=config['ordered_samples'])
    output:
        #bed = 'analysis/covtobed/co_intersect/covered.10x.merge.bed'
        bed = 'analysis/covtobed/co_intersect/covered.10x.target.bed'
    params:
        tempDir = 'analysis/covtobed/co_intersect/temp'
    shell:
        "python bin/parser.py ExtractCoIntersect"
        " --inputs {input.beds}"
        " --output {output.bed}"
        " --temp {params.tempDir}"

rule cnvkit_batch:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam',
        #bed = 'analysis/covtobed/{sample}/covered.10x.target.bed'
        #bed = 'analysis/covtobed/{sample}/covered.10x.merge.bed'
        #bed = 'analysis/covtobed/co_intersect/covered.10x.merge.bed'
        #bed = 'analysis/covtobed/co_intersect/covered.10x.target.bed'
    output:
        cnn = 'analysis/cnvkit/{sample}/{sample}.cnn',
        cnr = 'analysis/cnvkit/{sample}/{sample}.cnr',
        cns = 'analysis/cnvkit/{sample}/{sample}.cns'
    params:
        outdir = 'analysis/cnvkit/{sample}',
        genome = config['genome_fasta'],
        access = config['access_bed'],
        refflat = config['ref_flat'],
        target = config['target_bed']
    threads: 8
    log:
        stdout="analysis/cnvkit/{sample}/cnvkit_batch.log.stdout",
        stderr="analysis/cnvkit/{sample}/cnvkit_batch.log.stderr"
    shell:
        "cnvkit.py batch"
        " {input.bam}"
        " --normal"
        " --targets {params.target}"
        " --fasta {params.genome}"
        " --access {params.access}"
        " --annotate {params.refflat}"
        " --output-reference {output.cnn}"
        " --output-dir {params.outdir}"
        " -p {threads}"
        " > {log.stdout} 2> {log.stderr}"
        #" --targets {input.bed}"

#rule cnvkit_segment:
#    input:
#        cnr = 'analysis/cnvkit/{sample}/{sample}.cnr'
#    output:
#        cns = 'analysis/cnvkit/{sample}/{sample}.cns'
#    log:
#        stdout="analysis/cnvkit/{sample}/cnvkit_segment.log.stdout",
#        stderr="analysis/cnvkit/{sample}/cnvkit_segment.log.stderr"
#    shell:
#        "cnvkit.py segment"
#        " {input.cnr}"
#        " -o {output.cns}"
#        " > {log.stdout} 2> {log.stderr}"


rule cnvkit_call:
    input:
        cns = 'analysis/cnvkit/{sample}/{sample}.cns'
    output:
        cns = 'analysis/cnvkit/{sample}/{sample}.call.cns'
    params:
        m = 'threshold',
        t = '-1.1,-0.4,0.3,0.7'
    shell:
        "cnvkit.py call"
        " {input.cns}"
        " -m {params.m}"
        " -t={params.t}"
        " -o {output.cns}"

rule cnvkit_call_bed:
    input:
        cns = 'analysis/cnvkit/{sample}/{sample}.call.cns'
    output:
        bed = 'analysis/cnvkit/{sample}/{sample}.call.bed'
    shell:
        "cnvkit.py export bed"
        " {input.cns}"
        " --show all"
        " --label-genes"
        " -o {output.bed}"

rule cnvkit_call_anno:
    input:
        bed = 'analysis/cnvkit/{sample}/{sample}.call.bed'
    output:
        bed = 'analysis/cnvkit/{sample}/{sample}.call.anno.bed'
    params:
        bed = config['ipsc57_bed']
    shell:
        "bedtools intersect"
        " -wa -wb"
        " -a {input.bed}"
        " -b {params.bed}"
        " > {output.bed}"

rule cnvkit_gainloss:
    input:
        cnr = 'analysis/cnvkit/{sample}/{sample}.cnr',
        cns = 'analysis/cnvkit/{sample}/{sample}.cns'
    output:
        gainloss = 'analysis/cnvkit/{sample}/{sample}.gene.gainloss.txt'
    params:
        t = 0.4,
        m = 5
    shell:
        "cnvkit.py gainloss"
        " {input.cnr}"
        " -s {input.cns}"
        " -t {params.t}"
        " -m {params.m}"
        " > {output.gainloss}"

rule cnvkit_scatter:
    input:
        cnr = 'analysis/cnvkit/{sample}/{sample}.cnr',
        cns = 'analysis/cnvkit/{sample}/{sample}.cns'
    output:
        pdf = 'analysis/cnvkit/{sample}/{sample}-scatter.pdf',
        png = 'analysis/cnvkit/{sample}/{sample}-scatter.png',
    shell:
        "cnvkit.py scatter"
        " {input.cnr}"
        " -s {input.cns}"
        " -o {output.pdf} &&"
        "convert -density 150"
        " {output.pdf}"
        " -quality 90"
        " {output.png}"

rule cnvkit_diagram:
    input:
        cnr = 'analysis/cnvkit/{sample}/{sample}.cnr',
        cns = 'analysis/cnvkit/{sample}/{sample}.cns'
    output:
        pdf = 'analysis/cnvkit/{sample}/{sample}-diagram.pdf',
        png = 'analysis/cnvkit/{sample}/{sample}-diagram.png',
    shell:
        "cnvkit.py diagram"
        " {input.cnr}"
        " -s {input.cns}"
        " -o {output.pdf} &&"
        "convert -density 150"
        " {output.pdf}"
        " -quality 90"
        " {output.png}"

rule parse_cnvkit:
    input:
        callcns = expand('analysis/cnvkit/{sample}/{sample}.call.cns', sample=config['ordered_samples']),
        annobed = expand('analysis/cnvkit/{sample}/{sample}.call.anno.bed', sample=config['ordered_samples'])
    output:
        refbed = 'analysis/cnvkit/compbed/ref.ipsc57LBvsRef.bed',
        compbed = 'analysis/cnvkit/compbed/comp.ipsc57LBvsRef.bed'
    params:
        bed = config['ipsc57_bed'],
        samples = config['ordered_samples'],
        outdir = 'analysis/cnvkit/compbed'
    shell:
        "python bin/parser.py"
        " ParseCNVkit"
        " --samples {params.samples}"
        " --refbed {params.bed}"
        " --outdir {params.outdir}"
        
rule cnvkit_heatmap:
    input:
        callcns = expand('analysis/cnvkit/{sample}/{sample}.cns', sample=config['ordered_samples'])
    output:
        png = 'analysis/cnvkit/heatmap.png'
    shell:
        "cnvkit.py heatmap"
        " {input.callcns}"
        " -d"
        " -o {output.png}"

rule run_qualimap_bamqc:
    input:
        bam = 'analysis/inbam/{sample}/{sample}.bam',
        bai = 'analysis/inbam/{sample}/{sample}.bai'
    output:
        html = 'analysis/qualimap/{sample}/qualimapReport.html',
        outdir = directory('analysis/qualimap/{sample}')
    params:
        outdir = 'analysis/qualimap/{sample}',
        memsize = '24G',
        bed = config['target_bed']
    threads: 8
    shell:
        'qualimap bamqc'
        ' --java-mem-size={params.memsize}'
        ' -bam {input.bam}'
        #' -gff {params.bed}'
        ' -outdir {params.outdir}'
        ' -nt {threads}'

rule prep_qualimap_multibamqc:
    input:
        qm_path = expand('analysis/qualimap/{sample}', sample=config['ordered_samples'])
    output:
        txt = 'analysis/qualimap/multi/input_data.txt'
    run:
        import os
        outfh = open(output.txt, 'w')
        for qm_path in input.qm_path:
            qm_path = os.path.abspath(qm_path)
            sample_name = qm_path.split('/')[-1]
            group_name = qm_path.split('/')[-1].split('-')[0]
            items = [sample_name, qm_path, group_name]
            outfh.write("{0}\n".format('\t'.join(items)))
        outfh.close()

rule run_qualimap_multibamqc:
    input:
        txt = 'analysis/qualimap/multi/input_data.txt'
    output:
        html = 'analysis/qualimap/multi/multisampleBamQcReport.html'
    params:
        outdir = 'analysis/qualimap/multi',
        memsize = '24G',
        bed = config['target_bed']
    threads: 8
    shell:
        'qualimap multi-bamqc'
        ' --java-mem-size={params.memsize}'
        ' --data {input.txt}'
        ' -gff {params.bed}'
        ' -outdir {params.outdir}'








