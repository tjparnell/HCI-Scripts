# HCI-Scripts
Miscellaneous bioinformatics scripts written by me for myself and other people.

Mostly (all?) Perl scripts based on Bio::ToolBox.

## Annotation
Scripts that are about working with annotation, either as a file or from a 
database.

#### `delete_db_features.pl`
Delete features from a `Bio::DB::SeqFeature::Store` database

#### `dump_db_feature.pl`
Dump features from a `Bio::DB::SeqFeature::Store` (or other BioPerl) database 
as perl code using `Data::Dumper`.

#### `get_ucsc_transcript_lengths.pl`
Calculate transcript lengths from a UCSC style table

#### `gtf2gff3.pl`
Convert a GTF gene table to a GFF3 gene table

#### `print_chromosome_lengths.pl`
Prints out chromosome lengths from any database or indexed data file, 
such as bam, bigwig, bigbed, or useq file

#### `sum_repeat_family_counts.pl`
quick script to sum up the number of counts for repeat family members

## BamFile
Scripts for working with bam files.

#### `bin_bam_mapq.pl`
Make histogram table of all the mapping quality scores in a bam file.

#### `count_alignments.pl`
Count the number of alignments to each chromosome or reference sequence.

#### `pull_matching_alignments.pl`
Extracts reads from a bam file based on their name.

#### `sum_bam_isizes.pl`
Make a histogram table of all the reported insertion sizes for paired-end 
alignments.
 
## ChIPSeq
Scripts for working with ChIPSeq.

#### `ChIPNexus_bam_dedup.pl`
Remove single-end duplicates from a bam taking into account the ChIPNexus 
unique barcode present in the name.

#### `ChIPNexus_fastq_barcode_processer.pl`
Processes a ChIPNexus barcode fastq file by extracting the random barcode 
from the sequence read and appending it to the name. Required before 
alignment and de-duplication.

#### `delta_delta_calculator.pl`
Calculates a delta-delta-rpm signal from four data tables.

#### `summarize_narrowPeaks.pl`
Summarize scores and lengths from narrowPeak files.

 
## Clusters
Scripts for working with k-means cluster data files.

#### `print_clusters.pl`
Wrapper script for the command line Java TreeView program for quickly 
printing lots of heat maps out at once. Must provide your own .jtv if 
you do not want default colors.
 
## DataTables
Scripts for working with `Bio::ToolBox` data tables.

#### `create_summary_file.pl`
Quickly print out a summary file of a data table, much like what 
`get_relative_data` will do. Transposes the data table for easier 
graphing purposes.

#### `summarize_data_tables.pl`
Alternative script for summarizing data tables by collapsing all 
the rows into a single row. Does not transpose the table.

#### `transpose.pl`
Transposes (rotates 90 degrees) a data table. 

#### `two_table_calculator.pl`
Generates a third data table by mathematically combining the cells of 
two existing data tables. Crude matrix math.
 
## Nucleosome
Scripts for working with MNase Sequencing of nucleosomes.

#### `bam_count_MNase_ends.pl`
Count how many reads end in a typical MNase cut-site `[AT][AT]`. 
Requires a genome fasta to get the other base.

#### `run_nucleosome.pl`
Wrapper script to run Eran Segal's nucleosome prediction program.
 
## RNASeq
Scripts for working with RNASeq data

#### `parse_RNASeqMetrics.pl`
Parses one or more Picard's CollectRNASeqMetrics output files into a single 
comprehensive output table suitable for graphing.
 
## RT
Scripts for customizing Best Practical's Request Tracker software.

#### HCIWhiteListFilter.pm
#### email_bioinformaticshelp
#### index.html
#### requests
#### submission.html
 
## SomaticVariants
Scripts for working with somatic variant calling

#### `add_annot_table_freqs.pl`
Adds variant frequencies as an additional column to an annotated table 
exported using USeq VCFReporter.

#### `filter_vcf_by_FA.pl`
A script to filter somatic VCF files based on the FA (Frequency of Alternate) 
attribute. 

#### `get_transcript_lengths_for_variants.pl`
Custom script to merge in transcript and CDS lengths for Variants to assist in 
calculating mutations rates.

#### `pull_vcf_FA.pl`
Script to extract the fraction alternate (FA) values from somatic VCF files.

#### `summarize_exported_somaticVariant_table.pl`
Summarize Somatic variant tables into gene hit count tables.  
Assume tables are from USeq VCFReporter text output. 

#### `update_somaticVCF_attributes.pl`
A script to fix and standardize sample attributes in somatic VCF files. 
 
## WigFiles
Scripts for working with wig files.

#### `convert_FE_to_log2FE.pl`
Old script for converting a macs2 generated fold enrichment (FE) bedgraph 
file into log2 fold enrichment.

#### `file2big.pl`
Wrapper script to convert any appropriate text file into it's big equivalent, 
i.e wig to bigWig, bed to bigBed. 

#### `manipulate_wig.pl`
Multipurpose script to manipulate the values in any form of text wig file, 
including fixedStep, variableStep, or bedGraph. Multiple conversions, 
transformations, and statistics are available.

#### `summarize_coverage.pl`
Summarize fraction of coverage depth from 1 or more bedGraph files.

#### `trim_wig.pl`
Old script to cut low scores from a wig file.



