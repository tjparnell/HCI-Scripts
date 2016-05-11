# HCI-Scripts
Miscellaneous bioinformatics scripts written by me for myself and other people.

Some are fully documented, polished applications exported from a purging of Bio::ToolBox at version 1.40. These no longer fit the intent of the original 
package, but are too useful to discard.

Some are quickie one-off scripts thrown together for some specific need. Most 
illustrate the utility of using the Bio::ToolBox libraries, so if you need 
some inspiration....

Mostly (all?) are Perl scripts based on Bio::ToolBox, available at 
https://github.com/tjparnell/biotoolbox.


## Installation

These scripts do not need to be installed anywhere special. They can be run 
most anywhere. Most require certain Perl modules to be installed, for 
example the `Bio::ToolBox` libraries found at this Repository. As long as 
modules are installed in a PERL5 library available to your Perl installation 
(either the site directory or through `local::lib`), they should run. 
Scripts will not run if required modules are missing. You can install 
missing modules using your favorite CPAN utility.


## Organization

The scripts are organized into the following folders based loosely on their function. 

### Annotation
Scripts that are about working with genomic annotation, derived either from 
a file or from a database.

### BamFile
Scripts for working with bam files. Requires installation of `Bio::DB::Sam`.

### ChIPSeq
Scripts for working with ChIPSeq data.

### Clusters
Scripts for working with k-means cluster data files.

### DataTables
Scripts for working with tab-delimited text data tables using the 
`Bio::ToolBox` libraries.

### Nucleosome
Scripts for working with MNase Sequencing of nucleosomes.

### RT
Scripts for customizing Best Practical's Request Tracker software.

### Sequence
Scripts for working with genomic sequence.

### SomaticVariants
Scripts for working with somatic variant calling. Demonstrates use of 
working with VCF files.

### WigFiles
Scripts for working with and manipulating wig files.




