####################################
# Should be sourced by other scripts
####################################

export rawdata_dir=$PROJDIR/rawData
export procdata_dir=$PROJDIR/processedData
export results_dir=$PROJDIR/results

export input_dirname=d0_Input
export trimgalore_dirname=d1_TrimGalore
export fastqc_dirname=d2_FastQC
export star_dirname=d3_STAR
export htseq_dirname=d4_HTSeq
export gsea_dirname=d7_GSEAinput
export bam_dirname=d8_BAM

export input_dir=$procdata_dir/$input_dirname
export trimgalore_dir=$procdata_dir/$trimgalore_dirname
export fastqc_dir=$procdata_dir/$fastqc_dirname
export star_dir=$procdata_dir/$star_dirname
export htseq_dir=$procdata_dir/$htseq_dirname
export gsea_dir=$procdata_dir/$gsea_dirname
export bam_dir=$procdata_dir/$bam_dirname

# UPDATE: if original files were symlinked rather than uncompressed+concatenated,
# remove first pattern
export remove_patterns="$input_dirname/*.fastq  $trimgalore_dirname/*.fq  $star_dirname/*.sam"
