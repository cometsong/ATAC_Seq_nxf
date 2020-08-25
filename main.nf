#!/usr/bin/env nextflow

import Helpers
import Logos

logo = new Logo()
println logo.show()

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    // Configurable variable parameters specific to individual runs:
    params.fastqInputs = null // paths to input fastqs. Comma delimited (e.g., /path/to/R1,/path/to/R2).
    params.outputDir   = null // Base directory where output will be stored.
    params.bowtieIndex = null // absolute path to bowtie2 index.
    params.threads     = 1    // threads. (Must be numeric.)
    params.chain       = null // path to chain file OR null (e.g., /path/to/REF-to-CAST_EiJ.chain)
    params.tmpdir      = "$TMPDIR"

    //NOTE: See nextflow.config for ATACSupportFiles actual default!
    params.ATACSupportFiles = ''   // ATAC-seq tool images
}
setParamDefaults()

def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow -c path/to/params.cfg run ${workflow.projectDir}/${workflow.manifest.name} -profile sumner
            (The params.cfg file needs to have the following mandatory parameters
             OR they need to specified on the command line.)

    Mandatory:
        --fastqInputs           Paths to input fastqs. Comma delimited (e.g., /path/to/R1,/path/to/R2)
        --outputDir             Base directory where output will be stored.
        --bowtieIndex           Absolute path to Bowtie2 index file
        --threads               Threads (Must be numeric.)
        --chain                 Path to chain file (e.g., /path/to/REF-to-CAST_EiJ.chain)

    Optional:
        --ATACSupportFiles      ATAC-seq tool container images
        --tmpdir                Path to hold temporary files for software execution.
        --email                 The email address to send the pipeline report.
        -name                   Name for the pipeline run. If not specified Nextflow will
                                    automatically generate a random mnemonic.

        -profile                Environment config to use.
                                    [ choices: standard (local), sumner ]
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
////// Required parameters \\\\\\
if ( ! params.fastqInputs ) {
    exit 1, "Parameter ERROR: fastqInputs ($params.fastqInputs) must be a comma separated list of two _R1 and _R2 fastq filenames."
}
if ( ! params.outputDir ) {
    exit 1, "Parameter ERROR: Output directory parameter must be specified."
}
if ( ! file(params.outputDir).exists() ) {
    exit 1, "Parameter ERROR: Missing base root output directory ($params.outputDir) is not found: check if path is correct."
}

// Check bowtieIndex, append file extensions to check exists
if ( ! params.bowtieIndex ) {
    exit 1, "Parameter ERROR: bowtieIndex absolute path must be specified."
} else {
    bowtieFile = params.bowtieIndex + '.1.bt2'
    if ( ! file(bowtieFile).exists() ) {
        exit 1, "Parameter ERROR: bowtie2 index file for param ($params.bowtieIndex) does not exist."
    }
}

////// Optional params \\\\\\
if ( ! file(params.ATACSupportFiles).exists() ) {
    exit 1, "Parameter ERROR: ATACSupportFiles directory ($params.ATACSupportFiles) does not exist."
}
if ( params.chain || params.chain==null ) {
    chain = (params.chain==null || params.chain.length()==0) ? null : params.chain
    if (chain!=null) {
        chain = file(chain).exists() ? file(chain) : 'no-such-file!'
    }
}
else {
    exit 1, "Parameter ERROR: Chain file ($params.chain) does not exist."
}
// Check threads is numeric value
def number_vars = [ 'threads' ]
CheckParams.is_numeric(params, number_vars)


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
def summary = [:]
summary['Pipeline']         = workflow.manifest.name
summary['Description']      = workflow.manifest.description
if(workflow.revision) {
    summary['Pipeline Release'] = workflow.revision
}
summary['Run Name']         = workflow.runName
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
summary['Config Files']     = workflow.configFiles
summary['Command Line']     = workflow.commandLine
summary['Nextflow Info']    = "v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Workflow dir']     = workflow.projectDir
if(workflow.containerEngine) {
    summary['Container Engine'] = "$workflow.containerEngine"
    //summary['Containers'] = "$workflow.container" // only works with Docker. :(
}

// Pipeline Params:
summary['Parameters......']  = ''
if(params.email) {
    summary['.  Email']      = params.email
}
summary['.  Output dir']     = params.outputDir
summary['.  FASTQ inputs']   = params.fastqInputs
summary['.  Bowtie Index']   = params.bowtieIndex
summary['.  Threads']        = params.threads
summary['.  Chain']          = params.chain
summary['.  ATACSupportFiles'] = params.ATACSupportFiles
summary['.  TMP dir']        = params.tmpdir

summary['Run Start Time']    = workflow.start

//print summary header:
println Summary.show(summary)

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def timestamp = new Date().format("yyyyMMdd'T'hhmmSSS") // `date +"%Y%M%dT%H%M%N"`
def tmpdir    = params.tmpdir

//~~~~~ Final Processed File Extension ~~~~~
def processed_bam_extension = "sorted.rmDup.rmChrM.rmMulti.filtered.shifted.mm10" //NOTE: add .bam or .ba*
//%s;\.sorted\.rmDup\.rmChrM\.rmMulti\.filtered\.shifted\.mm10;.${processed_bam_extension};g

//~~~~~ Determine sampleID and outputDir ~~~~
//FIXME: Here we presume '_R1' is used in fastq name! (e.g. not '-R1', '.R1', etc)
def fqPairRE = ~/_R1.*\..*f[ast]*q.*$/
fqin = params.fastqInputs.tokenize(",")
fqR1 = file(fqin[0])

fqR1p = fqR1.toAbsolutePath().toString() // bc relative paths cause symlink breakage
fqR2p = file(fqin[1]).toAbsolutePath().toString()

def sampleID = ( fqR1.name - fqPairRE )

//~~~~~~~~~~ Initial Channel of SampleID and Fastq Data ~~~~
Channel.of( sampleID, fqR1p, fqR2p )
       .toList()
       .set { sample_fastqs_ch }

//~~~~~~~~~~~~~~~~ Primary publishDir == sample_outdir ~~~~~
def sample_outdir = "${params.outputDir}/${timestamp}_${sampleID}/"


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Processes ~ ~ ~ ~ ~
/*
Trim FASTQ files	(cutadapt)
QC of FASTQ	(fastqc)
Align FASTQ files	(bowtie)
Sort the alignment	(samtools)
Flagging PCR Duplicate	(gatk)
Mark duplicates, remove them	(samtools)
Calculate %mtDNA and Filter Mitochondrial Reads	(samtools)
Filter Non-Unique and Include Only 'properly mapped reads' Alignments	(samtools)
Shifting Reads	(samtools)
    "alignmentSieve"	(deeptools)
    Re-sort BAM alignment	(samtools)
IF (chain):
    Convert Peak Coordinates to B6	(g2gtools)
    Sort BAM by coordinates	(samtools)
    Validate SAM File	(gatk)
    Filter Bad SAM Reads	(gatk)
    Re-sort BAM by name	(samtools)
    Fix Mate TLEN Info	(samtools)
    Re-sort BAM by coordinates	(samtools)
    	(samtools: index, view, grep, index)
ELSE (no chain):
   	(samtools: index, view, grep, index)
Peak Calling	(macs2)
Fraction of reads in peaks FRiP	(samtools, bedtools, awk)
Get coverage in each peak	(awk)
Feature counts	(subread)
Quality checks:
    Fragment/Insert size	(samtools, fragment_length_plot.R)
    Library Complexity	(samtools)
    Calculate PBC metrics	(bedtools)
LogParser
*/


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Trim Cutadapt ~~~~~
// Check parameters are numbers
number_vars = [ 'cutadaptMinLength', 'cutadaptQualCutoff', ]
CheckParams.is_numeric(params, number_vars)
process trim_fastq {
    tag "$sampleID"
    label 'cutadapt'
    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy' // also to LogParser
    cpus params.threads

    input:
    tuple sampleID, fqR1, fqR2 from sample_fastqs_ch
    output:
    tuple sampleID, file("${sampleID}_R{1,2}_paired_trimmed.fq") \
          into ( trimmed_reads_fq, trim_fastqc )
    tuple sampleID, file("*.log") into log_cutadapt

    script:
    log.info "-----Cutadapt running on ${sampleID}-----"
    """
    cutadapt \
        -a ${params.cutadaptAdapterR1} \
        -A ${params.cutadaptAdapterR2} \
        --minimum-length ${params.cutadaptMinLength} \
        --quality-cutoff ${params.cutadaptQualCutoff} \
        -j ${params.threads} \
        -o ${sampleID}_R1_paired_trimmed.fq \
        -p ${sampleID}_R2_paired_trimmed.fq \
        ${fqR1} ${fqR2} \
        > ${sampleID}_cutadapt.log \
        2>&1
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FastQC ~~~~~
process fastqc {
    tag "$sampleID"
    label 'fastqc'
    publishDir "${sample_outdir}", pattern: "*_fastqc.*", mode: 'move'
    cpus params.threads

    input:
    tuple sampleID, file(trimmedfqs) from trim_fastqc
    output:
    file("*_fastqc.*")

    script:
    log.info "-----FastQC running on ${sampleID}-----"
    """
    fastqc -t ${params.threads} ${trimmedfqs[0]} ${trimmedfqs[1]}
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Align Trimmed Bowtie ~~~~~
// Check parameters are numbers
number_vars = [ 'bowtieMaxInsert' ]
CheckParams.is_numeric(params, number_vars)
process align_trimmed_fastq {
    tag "$sampleID"
    label 'bowtie2'
    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy' // also to LogParser
    cpus params.threads

    input:
    tuple sampleID, file(trimmedfqs) from trimmed_reads_fq
    output:
    tuple sampleID, file("${sampleID}.sam") into aligned_ch
    file("${sampleID}_bowtie2.log") into log_bowtie

    script:
    log.info "-----Bowtie running on ${sampleID}-----"
    """
    bowtie2 \
        --very-sensitive \
        -X ${params.bowtieMaxInsert} \
        -q \
        -p ${params.threads} \
        -x ${params.bowtieIndex} \
        -1 ${sampleID}_R1_paired_trimmed.fq \
        -2 ${sampleID}_R2_paired_trimmed.fq \
        -S ${sampleID}.sam \
        2>${sampleID}_bowtie2.log
    """
}
process sort_alignment {
    tag "$sampleID"
    label 'samtools'

    input:
    tuple sampleID, file(sample_aligned) from aligned_ch
    output:
    tuple sampleID, file("${sampleID}.sorted.bam*") into sorted_align_ch

    script:
    log.info "-----Sorting and indexing Bowtie alignment of ${sampleID}-----"
    """
    samtools sort \
        -@ ${params.threads} \
        -O bam \
        -o ${sampleID}.sorted.bam \
        ${sample_aligned}
    samtools index \
        ${sampleID}.sorted.bam
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Flagging PCR Duplicates ~~~~~
process flag_pcr_dupes {
    tag "$sampleID"
    label 'gatk'
    publishDir "${sample_outdir}", pattern: "*.sorted.metrics", mode: 'copy' // also to LogParser
    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy' // also to LogParser
    cpus params.threads

    input:
    tuple sampleID, file(sorted_bams) from sorted_align_ch
    output:
    tuple sampleID, file("${sampleID}.sorted.marked.ba*") \
          into ( sorted_marked_bam_rmdup_ch, sorted_marked_bam_lib_cmplx_ch )
    file("${sampleID}.sorted.metrics") into log_sorted_metrics
    file("${sampleID}.picard.log") into log_picard

    script:
    log.info "-----Flagging PCR Duplicates on ${sampleID}-----"
    """
    gatk \
        --java-options "-XX:ParallelGCThreads=${params.threads} -Djava.io.tmpdir=${tmpdir}" \
        MarkDuplicates \
        --INPUT ${sorted_bams[0]} \
        --OUTPUT ${sampleID}.sorted.marked.bam \
        --METRICS_FILE ${sampleID}.sorted.metrics \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR ${tmpdir} \
        > ${sampleID}.picard.log 2>&1
    """
}
process rm_dupe_reads {
    tag "$sampleID"
    label 'samtools'

    input:
    tuple sampleID, file(marked_bams) from sorted_marked_bam_rmdup_ch
    output:
    tuple sampleID, file("${sampleID}.sorted.rmDup.ba*") into sorted_rmdup_bam_ch

    script:
    """
    samtools view -h -b -F 1024 ${marked_bams[1]} > ${sampleID}.sorted.rmDup.bam
    samtools index ${sampleID}.sorted.rmDup.bam
    """
}


//~~~ Calculate %mtDNA and Filter Mitochondrial Reads ~~~~~
process calc_mtdna_filter_chrm {
    tag "$sampleID"
    label 'samtools'
    publishDir "${sample_outdir}", pattern: "*_mtDNA_Content.txt", mode: 'copy' // also to LogParser
    cpus params.threads

    input:
    tuple sampleID, file(rmdup_bams) from sorted_rmdup_bam_ch
    output:
    tuple sampleID, file("${sampleID}.sorted.rmDup.rmChrM.ba*") into mtdna_filter_bams_ch
    file("${sampleID}_mtDNA_Content.txt") into log_mtdna_content

    shell:
    log.info "-----Calculate %mtDNA and Filter Mitochondrial Reads on ${sampleID}-----"
    '''
    mtReads=$(samtools idxstats !{rmdup_bams[0]} | grep 'MT' | cut -f 3)
    totalReads=$(samtools idxstats !{rmdup_bams[0]} | awk '{SUM += $3} END {print SUM}')

    echo 'mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > !{sampleID}_mtDNA_Content.txt

    samtools view -@ !{params.threads} -h !{rmdup_bams[0]} \
    | grep -v MT \
    | samtools sort -@ !{params.threads} -O bam \
        -o !{sampleID}.sorted.rmDup.rmChrM.bam \
    && samtools index !{sampleID}.sorted.rmDup.rmChrM.bam
    '''
}


//~~~~~Filter Non-Unique and Include Only 'properly mapped reads' Alignments~~~~~
process filter_rmmulti_shift {
    tag "$sampleID"
    label 'samtools'
    publishDir "${sample_outdir}", pattern: "*.sorted.rmDup.rmChrM.rmMulti.filtered.ba*", mode: 'copy'
    cpus params.threads

    input:
    tuple sampleID, file(mtdna_bams) from mtdna_filter_bams_ch
    output:
    tuple sampleID, file("${sampleID}.shift.tmp0.ba*") into shift_tmp_bams_ch
    tuple sampleID, file("${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.ba*") \
          into sort_rm_dup_chrm_multi_filter_ch

    script:
    log.info "-----Filter Non-Unique and Include Only 'properly mapped reads' Alignments on ${sampleID}-----"
    """
    # Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
    # Retain properly paired reads -f 2
    samtools view \
        -@ ${params.threads} -h -q 30 ${mtdna_bams[0]} \
        > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.bam
    samtools view \
        -@ ${params.threads} -h -b -F 1804 -f 2 \
        ${sampleID}.sorted.rmDup.rmChrM.rmMulti.bam \
        > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam
    samtools index \
        ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam
    samtools sort \
        -@ ${params.threads} -O bam \
        -o ${sampleID}.shift.tmp0.bam \
        ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam
    samtools index \
        ${sampleID}.shift.tmp0.bam
    """
}
process filter_rmmulti_sieve {
    tag "$sampleID"
    label 'deeptools'
    cpus params.threads

    input:
    tuple sampleID, file(shift_bams) from shift_tmp_bams_ch
    output:
    tuple sampleID, file("${sampleID}.shift.tmp.ba*") into shift_bams_ch

    script:
    log.info "-----Running deeptools alignmentSieve on ${sampleID}-----"
    """
    alignmentSieve \
        --numberOfProcessors ${params.threads} \
        --ATACshift \
        --bam ${shift_bams[0]} \
        -o ${sampleID}.shift.tmp.bam
    """
}
process filter_rmmulti_sort {
    tag "$sampleID"
    label 'samtools'
    cpus params.threads

    input:
    tuple sampleID, file(shift_bams) from shift_bams_ch
    output:
    tuple sampleID, file("${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.ba*") \
        into ( sort_rm_dup_chrm_mult_filt_shift_chain_file,
               sort_rm_dup_chrm_mult_filt_shift_chain_null )

    script:
    log.info "-----Re-sorting shifted bam on ${sampleID}-----"
    """
    samtools sort \
        -@ ${params.threads} -O bam \
        -o ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.bam \
        ${shift_bams[0]}
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ IF chain ~~~~~
//~~~~~~~~~~~~~~~~~~~~~ Convert Peak Coordinates to B6 ~~~~~
process chain_convert_peak_b6 {
    tag "$sampleID"
    label 'g2gtools'

    input:
    tuple sampleID, file(bam_shifted) from sort_rm_dup_chrm_mult_filt_shift_chain_file
    output:
    tuple sampleID, file("${sampleID}.tmp.mm10.ba*") into bam_peak_b6_mm10_ch

    when: chain != null

    script:
    log.info "-----Convert Peak Coordinates to B6 on ${sampleID}-----"
    """
    # convert to mm10 coordinates
    g2gtools convert \
        -r -f bam -c ${chain} \
        -i ${bam_shifted[0]} \
        -o ${sampleID}.tmp.mm10.bam
    """
}
process chain_sort_coords {
    tag "$sampleID"
    label 'samtools'
    cpus params.threads

    input:
    tuple sampleID, file(bam_peak_b6_mm10) from bam_peak_b6_mm10_ch
    output:
    tuple sampleID, file("${sampleID}.tmp1.mm10.ba*") \
        into ( bams_sort_mm10_extract_ch, bams_sort_mm10_filter_ch )

    when: chain != null

    script:
    """
    # sort bam by coordinates
    samtools sort \
        -@ ${params.threads} -O bam \
        -o ${sampleID}.tmp1.mm10.bam \
        ${bam_peak_b6_mm10[0]}
    """
}
process chain_extract_badreads_b6 {
    tag "$sampleID"
    label 'gatk'
    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'
    validExitStatus 0,3,4 //,o,t,h,e,r,s?

    input:
    tuple sampleID, file(bam_sort_mm10) from bams_sort_mm10_extract_ch
    output:
    tuple sampleID, file("BAD_READS") into bad_reads_ch
    file("${sampleID}_ValidateSamFile.log") 

    when: chain != null

    script:
    """
    # extract a list of 'bad reads'
    gatk ValidateSamFile \
        -I ${bam_sort_mm10[0]} \
        -MODE VERBOSE -MO 10000000 \
        -O BAD_READS \
        --IGNORE MISSING_READ_GROUP \
        --IGNORE RECORD_MISSING_READ_GROUP \
        --IGNORE MISSING_TAG_NM \
        --IGNORE QUALITY_NOT_STORED \
        > ${sampleID}_ValidateSamFile.log 2>&1
    """
}
process chain_bad_to_uniq_reads {
    tag "$sampleID"
    label 'utils'

    input:
    tuple sampleID, file(bad_reads) from bad_reads_ch
    output:
    file("ReadName_unique") into read_uniq_ch

    when: chain != null

    shell:
    '''
    # Remove 'bad reads' from bam file
    cat !{bad_reads} \
    | awk '{print $5}' \
    | sed -r 's/\\,//g' \
    | sort -n \
    | uniq -c \
    | awk '{print $2}' \
    > ReadName_unique
    '''
}
process chain_filter_reads_b6 {
    tag "$sampleID"
    label 'gatk'
    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'

    input:
    tuple sampleID, file(bam_sort_mm10) from bams_sort_mm10_filter_ch
    file(ReadName_unique) from read_uniq_ch
    output:
    tuple sampleID, file("${sampleID}.tmp2.mm10.ba*") into bam_mm10_ch
    file("${sampleID}_FilterSamReads.log")

    when: chain != null

    script:
    """
    # filter list to unique names (reads can fail for multiple reasons appear on list more than once)
    gatk FilterSamReads \
        -I ${bam_sort_mm10[0]} \
        -RLF ReadName_unique \
        --FILTER excludeReadList \
        --VALIDATION_STRINGENCY LENIENT \
        -O ${sampleID}.tmp2.mm10.bam \
        > ${sampleID}_FilterSamReads.log 2>&1
    """
}
process chain_sort_fixmate_bam {
    tag "$sampleID"
    label 'samtools'
    publishDir "${sample_outdir}", pattern: "*.${processed_bam_extension}.ba*", mode: 'copy'
    cpus params.threads

    input:
    tuple sampleID, file(bam_mm10) from bam_mm10_ch
    output:
    tuple sampleID, file("${sampleID}.${processed_bam_extension}.ba*") \
          into sort_rm_dup_chrm_multi_filter_shifted_mm10_chain

    when: chain != null

    script:
    """
    # name sort the bam
    samtools sort \
        -n \
        -@ ${params.threads} -O bam \
        -o ${sampleID}.tmp3.mm10.bam ${bam_mm10[0]}
    # fix the mate information. This is done to fix 'TLEN' which is required for MACS2.
    samtools fixmate \
        -O bam ${sampleID}.tmp3.mm10.bam ${sampleID}.tmp4.mm10.bam
    # re-sort bam by coordinates.
    samtools sort \
        -@ ${params.threads} -O bam \
        -o ${sampleID}.tmp5.mm10.bam ${sampleID}.tmp4.mm10.bam
    samtools index ${sampleID}.tmp5.mm10.bam
    samtools view ${sampleID}.tmp5.mm10.bam \
        -h 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
        | grep -ve 'SN:MT*' > tmp.sam
    samtools view -b tmp.sam \
        > ${sampleID}.${processed_bam_extension}.bam
    ### Note, for some reason the command that was working for mm10 only samples as seen
    ### below was not working here. This work around provided a valid BAM file
    samtools index \
        ${sampleID}.${processed_bam_extension}.bam
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ NOT chain ~~~~~
process non_chain_reindex {
    tag "$sampleID"
    label 'samtools'
    publishDir "${sample_outdir}", pattern: "*.${processed_bam_extension}.ba*", mode: 'copy'

    input:
    tuple sampleID, file(bam_shifted) from sort_rm_dup_chrm_mult_filt_shift_chain_null
    output:
    tuple sampleID, file("${sampleID}.${processed_bam_extension}.ba*") \
          into sort_rm_dup_chrm_multi_filter_shifted_mm10_non_chain

    when: chain == null

    script:
    """
    samtools index ${bam_shifted[0]}
    samtools view  ${bam_shifted[0]} \
        -h 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
        | grep -ve 'SN:MT*\\|SN:GL*\\|SN:JH:*' > tmp.sam
    samtools view -b tmp.sam \
      > ${sampleID}.${processed_bam_extension}.bam
    # re-header to remove unplaced/unloacalized
    samtools index \
      ${sampleID}.${processed_bam_extension}.bam
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Final Processing ~~~~~

//~~~~~~~~~~~~~~~~~ Channel Check [Non|Chain] Complete ~~~~~
sort_rm_dup_chrm_multi_filter_shifted_mm10_chain
.mix( sort_rm_dup_chrm_multi_filter_shifted_mm10_non_chain )
.into{
    processed_bams_ch_peak_calling;
    processed_bams_ch_bigwig;
    processed_bams_ch_frip_reads;
    processed_bams_ch_final_calc_frip;
    processed_bams_ch_feature_counts;
    }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Peak Calling ~~~~~
process peak_calling {
    tag "$sampleID"
    label 'macs2'
    publishDir "${sample_outdir}", pattern: "*_peaks.narrowPeak", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*_summits.bed", mode: 'move'
    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'

    scratch "${tmpdir}"

    input:
    tuple sampleID, file(processed_bams) from processed_bams_ch_peak_calling
    output:
    tuple sampleID, file("${sampleID}_peaks.narrowPeak") \
        into ( peak_frip_reads_ch, peak_call_cvg_ch )
    file("${sampleID}_summits.bed")
    file("${sampleID}_macs2.log")

    script:
    log.info "-----Peak Calling on ${sampleID}-----"
    """
    # -f BAMPE, use paired-end information
    # --keep-dup all, keep all duplicate reads.
    macs2 callpeak \
        -f BAMPE \
        --nomodel \
        -g mm \
        --keep-dup all \
        --cutoff-analysis \
        --tempdir ${tmpdir} \
        -n ${sampleID} \
        -t ${processed_bams[0]} \
        --outdir . \
        > ${sampleID}_macs2.log 2>&1
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Peak Calling ~~~~~
process bam_coverage_bigwig {
    tag "$sampleID"
    label 'deeptools'
    publishDir "${sample_outdir}", pattern: "*.bigwig", mode: 'copy'
    cpus params.threads

    input:
    tuple sampleID, file(processed_bams) from processed_bams_ch_bigwig
    output:
    file("*.bigwig")

    script:
    log.info "-----Running deeptools bamCoverage bigwig on ${sampleID}-----"
    """
    bamCoverage \
        --numberOfProcessors ${params.threads} \
        --binSize 10 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2652783500 \
        --bam *filtered.shifted.mm10.bam \
        --outFileFormat bigwig \
        --outFileName ${sampleID}.bigwig
    """
}


//~~~~~~~~~~~~~~~~~~ Fraction of reads in peaks (FRiP) ~~~~~
process frip_reads_in_peaks {
    tag "$sampleID"
    label 'bedtools'

    input:
    tuple sampleID, file(narrow_peaks) from peak_frip_reads_ch
    tuple sampleID, file(processed_bams) from processed_bams_ch_frip_reads
    output:
    tuple sampleID, file("reads_in_peaks.tmp.ba*") into reads_peaks_bams_ch

    script:
    log.info "-----Fraction of reads in peaks (FRiP) on ${sampleID}-----"
    """
    bedtools sort \
        -i ${narrow_peaks} \
    | bedtools merge -i stdin \
    | bedtools intersect -u -nonamecheck \
        -a ${processed_bams[0]} \
        -b stdin \
        -ubam \
    > reads_in_peaks.tmp.bam
    """
}
process final_calc_frip {
    tag "$sampleID"
    label 'samtools'
    publishDir "${sample_outdir}", pattern: "*_Fraction_reads_in_peak.txt", mode: 'copy' // also to LogParser

    input:
    tuple sampleID, file(reads_peaks_bams) from reads_peaks_bams_ch
    tuple sampleID, file(processed_bams)   from processed_bams_ch_final_calc_frip
    output:
    file("${sampleID}_Fraction_reads_in_peak.txt") into log_fraction_reads

    shell:
    '''
       total_reads=$(samtools view -c !{processed_bams[0]})
    reads_in_peaks=$(samtools view -c !{reads_peaks_bams[0]})
    FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
    echo -e ${FRiP}"\\t"${total_reads} \
        > !{sampleID}_Fraction_reads_in_peak.txt
    '''
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~ Get coverage in each peak ~~~~~
process peak_coverage {
    tag "$sampleID"
    label 'utils'
    publishDir "${sample_outdir}", pattern: "*_peaks.narrowPeak.*", mode: 'copy'

    input:
    tuple sampleID, file(narrow_peaks) from peak_call_cvg_ch
    output:
    tuple sampleID, file("${sampleID}_peaks.narrowPeak.saf") into peak_cvg_saf_ch

    shell:
    log.info "-----Get coverage in each peak on ${sampleID}-----"
    '''
    ## Make SAF file (+1 because SAF is 1-based, BED/narrowPeak is 0-based)
    awk 'OFS="\\t" {print $1"."$2"."$3, $1, $2, $3, "."}' !{narrow_peaks} \
    > !{sampleID}_peaks.narrowPeak.saf
    '''
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Feature Counts ~~~~~
process feature_counts {
    tag "$sampleID"
    label 'subread'
    cpus params.threads
    publishDir "${sample_outdir}", pattern: "*_peaks_countMatrix.txt", mode: 'copy'

    input:
    tuple sampleID, file(peak_cvg_saf) from peak_cvg_saf_ch
    tuple sampleID, file(processed_bams) from processed_bams_ch_feature_counts
    output:
    tuple sampleID, file("${sampleID}_peaks_countMatrix.txt") into feat_count_ch

    script:
    log.info "-----Feature Counts on ${sampleID}-----"
    """
    featureCounts \
        -a ${peak_cvg_saf} \
        -F SAF -p \
        -T ${params.threads} \
        -o ${sampleID}_peaks_countMatrix.txt \
        ${processed_bams[0]}
    """
}
process feature_count_to_bed {
    tag "$sampleID"
    label 'utils'
    publishDir "${sample_outdir}", pattern: "*_peaks_countMatrix.mm10.bed", mode: 'copy'

    input:
    tuple sampleID, file(peak_cnt_matrx) from feat_count_ch
    output:
    file("${sampleID}_peaks_countMatrix.mm10.bed")

    shell:
    '''
    tail -n +3 !{peak_cnt_matrx} \
        | awk -F $'\\t' 'BEGIN {OFS = FS} { print $2, $3, $4, $7, $6 }' \
        > !{sampleID}_peaks_countMatrix.mm10.bed
    '''
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Quality checks ~~~~~
process quality_checks {
    tag "$sampleID"
    label 'samtools'
    publishDir "${sample_outdir}", pattern: "*.fragment_length_count.txt", mode: 'copy'
    cpus params.threads

    input:
    tuple sampleID, file(sort_rm_filter_bam) from sort_rm_dup_chrm_multi_filter_ch
    output:
    tuple sampleID, file("${sampleID}.fragment_length_count.txt") into frag_len_count_ch

    script:
    log.info "-----Quality checks on ${sampleID}-----"
    log.info "-----Fragment/Insert size on ${sampleID}-----"
    """
    samtools view \
        -@ ${params.threads} ${sort_rm_filter_bam[0]} \
        | awk '\$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n \
        | sed -e 's/^[ \\t]*//' > ${sampleID}.fragment_length_count.txt
    """
}
process frag_len_plot {
    tag "$sampleID"
    label 'r_lang'
    publishDir "${sample_outdir}", pattern: "*fraglen_plot.pdf", mode: 'move', \
        saveAs: { pdf -> "${sampleID}.${pdf}" }

    input:
    tuple sampleID, file(frag_len_count) from frag_len_count_ch
    output:
    file("*fraglen_plot.pdf")

    script:
    log.info "-----Fragment Length Plot on ${sampleID}-----"
    """
    Rscript ${projectDir}/bin/fragment_length_plot.R ${frag_len_count}
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Library Complexity ~~~~~
process library_complexity {
    tag "$sampleID"
    label 'samtools'
    cpus params.threads

    input:
    tuple sampleID, file(sorted_marked_bam) from sorted_marked_bam_lib_cmplx_ch
    output:
    tuple sampleID, file("tmp.ba*") into tmp_bams_ch

    script:
    log.info "-----Library Complexity on ${sampleID}-----"
    //NOTE: the bam used here is sorted bam after duplicates marking sort bam by names
    """
    samtools sort \
        -@ ${params.threads} \
        -n \
        -O BAM \
        -o tmp.bam \
        ${sorted_marked_bam[1]}
    """
}
process calc_pbc_metrics {
    tag "$sampleID"
    label 'bedtools'
    publishDir "${sample_outdir}", pattern: "*.pbc.qc", mode: 'copy' // also to LogParser

    input:
    tuple sampleID, file(tmp_bams) from tmp_bams_ch
    output:
    file("${sampleID}.pbc.qc") into log_pbc_qc

    shell:
    log.info "-----Calculate PBC Metrics on ${sampleID}-----"
    '''
    bedtools bamtobed -bedpe \
        -i !{tmp_bams[0]} \
        | awk 'BEGIN{OFS="\\t"}{print $1,$2,$4,$6,$9,$10}' \
        | grep -v 'MT' | sort | uniq -c \
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} \
        ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
        END{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
        > !{sampleID}.pbc.qc
    '''
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LogParser ~~~~~
process logparser {
    tag "$sampleID"
    label 'logparser'
    publishDir "${sample_outdir}", pattern: "*.summary_QC_metrics.txt", mode: 'move'

    input:
    tuple sampleID, file(log_cutadapt) from log_cutadapt
                    file(log_bowtie)
                    file(log_sorted_metrics)
                    file(log_mtdna_content)
                    file(log_pbc_qc)
                    file(log_fraction_reads)
    output:
    file("${sampleID}.summary_QC_metrics.txt")

    script:
    log.info "-----LogParser on ${sampleID}-----"
    """
    python /LogParser.py > ${sampleID}.summary_QC_metrics.txt
    """
}


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Summary ~ ~ ~ ~ ~
workflow.onComplete {
    wfEnd = [:]
    wfEnd['Completed at'] = workflow.complete
    wfEnd['Duration']     = workflow.duration
    wfEnd['Exit status']  = workflow.exitStatus
    wfEnd['Success']      = workflow.success
    if(!workflow.success){
        wfEnd['!!Execution failed'] = ''
        wfEnd['.    Error']   = workflow.errorMessage
        wfEnd['.    Report']  = workflow.errorReport
    }
    Summary.show(wfEnd)
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// vim: set ft=groovy.nextflow ts=4 sw=0 tw=100 et fdm=syntax:
