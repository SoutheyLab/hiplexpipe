'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from, regex
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='hiplexpipe')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Call variants using undr_rover
    pipeline.transform(
        task_func=stages.apply_undr_rover,
        name='apply_undr_rover',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq.gz'),
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_R2_{lib[0]}.fastq.gz'),
        extras=['{sample[0]}'],
        output='variants/undr_rover/{sample[0]}.vcf')

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_R2_{lib[0]}.fastq.gz'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}', '{lib[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}.clipped.bam')

    # Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.clipped.bam'),
        output='.clipped.sort.bam')

    # High quality and primary alignments
    pipeline.transform(
        task_func=stages.primary_bam,
        name='primary_bam',
        input=output_from('sort_bam_picard'),
        filter=suffix('.clipped.sort.bam'),
        output='.clipped.sort.hq.bam')

    # index bam file
    pipeline.transform(
        task_func=stages.index_sort_bam_picard,
        name='index_bam',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.clipped.sort.hq.bam.bai')

    # generate mapping metrics.
    pipeline.transform(
        task_func=stages.generate_amplicon_metrics,
        name='generate_amplicon_metrics',
        input=output_from('primary_bam'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).clipped.sort.hq.bam'),
        output='alignments/metrics/{sample[0]}.amplicon-metrics.txt',
        extras=['{sample[0]}'])

    pipeline.transform(
        task_func=stages.intersect_bed,
        name='intersect_bed',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.intersectbed.bam')

    pipeline.transform(
        task_func=stages.coverage_bed,
        name='coverage_bed',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.bedtools_hist_all.txt')

    pipeline.transform(
        task_func=stages.genome_reads,
        name='genome_reads',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.mapped_to_genome.txt')

    pipeline.transform(
        task_func=stages.target_reads,
        name='target_reads',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.mapped_to_target.txt')

    pipeline.transform(
        task_func=stages.total_reads,
        name='total_reads',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.total_raw_reads.txt')

    pipeline.collate(
        task_func=stages.generate_stats,
        name='generate_stats',
        input=output_from('coverage_bed', 'genome_reads', 'target_reads', 'total_reads'),
        filter=regex(r'.+/(.+BS\d{4,6}.+S\d+)\..+\.txt'),
        output=r'all_sample.summary.\1.txt',
        extras=[r'\1', 'all_sample.summary.txt'])

    ###### GATK VARIANT CALLING ######
    # Call variants using GATK
    (pipeline.transform(
        task_func=stages.call_haplotypecaller_gatk,
        name='call_haplotypecaller_gatk',
        input=output_from('primary_bam'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9-_]+).clipped.sort.hq.bam'),
        output='variants/gatk/{sample[0]}.g.vcf')
        .follows('index_sort_bam_picard'))

    # Combine G.VCF files for all samples using GATK
    pipeline.merge(
        task_func=stages.combine_gvcf_gatk,
        name='combine_gvcf_gatk',
        input=output_from('call_haplotypecaller_gatk'),
        output='variants/gatk/ALL.combined.vcf')

    # Genotype G.VCF files using GATK
    pipeline.transform(
        task_func=stages.genotype_gvcf_gatk,
        name='genotype_gvcf_gatk',
        input=output_from('combine_gvcf_gatk'),
        filter=suffix('.combined.vcf'),
        output='.raw.vcf')

    # Apply GT filters to genotyped vcf
    pipeline.transform(
        task_func=stages.genotype_filter_gatk,
        name='genotype_filter_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.raw.vcf'),
        output='.raw.gt-filter.vcf')

    # Decompose and normalise multiallelic sites
    pipeline.transform(
        task_func=stages.vt_decompose_normalise,
        name='vt_decompose_normalise',
        input=output_from('genotype_filter_gatk'),
        filter=suffix('.raw.gt-filter.vcf'),
        output='.raw.gt-filter.decomp.norm.vcf')

    # Annotate VCF file using GATK
    pipeline.transform(
        task_func=stages.variant_annotator_gatk,
        name='variant_annotator_gatk',
        input=output_from('vt_decompose_normalise'),
        filter=suffix('.raw.gt-filter.decomp.norm.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.vcf')

    # Filter vcf
    pipeline.transform(
        task_func=stages.gatk_filter,
        name='gatk_filter',
        input=output_from('variant_annotator_gatk'),
        filter=suffix('.raw.gt-filter.decomp.norm.annotate.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.filter.vcf')

     #Apply VEP
    (pipeline.transform(
        task_func=stages.apply_vep,
        name='apply_vep',
        input=output_from('gatk_filter'),
        filter=suffix('.raw.gt-filter.decomp.norm.annotate.filter.vcf'),
        add_inputs=add_inputs(['variants/undr_rover/combined_undr_rover.vcf.gz']),
        output='.raw.gt-filter.decomp.norm.annotate.filter.vep.vcf')
        .follows('index_final_vcf'))

    #### concatenate undr_rover vcfs ####
    pipeline.transform(
        task_func=stages.sort_vcfs,
        name='sort_vcfs',
        input=output_from('apply_undr_rover'),
        filter=suffix('.vcf'),
        output='.sorted.vcf.gz')

    pipeline.transform(
        task_func=stages.index_vcfs,
        name='index_vcfs',
        input=output_from('sort_vcfs'),
        filter=suffix('.sorted.vcf.gz'),
        output='.sorted.vcf.gz.tbi')

    (pipeline.merge(
        task_func=stages.concatenate_vcfs,
        name='concatenate_vcfs',
        input=output_from('sort_vcfs'),
        output='variants/undr_rover/combined_undr_rover.vcf.gz')
        .follows('index_vcfs'))

    pipeline.transform(
        task_func=stages.index_final_vcf,
        name='index_final_vcf',
        input=output_from('concatenate_vcfs'),
        filter=suffix('.vcf.gz'),
        output='.vcf.gz.tbi')

    return pipeline
