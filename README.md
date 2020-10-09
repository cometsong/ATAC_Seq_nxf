## ATAC_Seq analysis pipeline (in Nextflow)

The Assay for Transposase-Accessible Chromatin followed by sequencing (ATAC-seq)
experiment provides genome-wide profiles of chromatin accessibility. Briefly,
the ATAC-seq method works as follows: loaded Tn5 transposase inserts sequencing
primers into open chromatin sites across the genome, and reads are then
sequenced. The ends of the reads mark open chromatin sites.

The developed pipeline maps trimmed paired-end Illumina reads from mouse
strains, to strain specific references. Mapped reads are then filtered to remove
duplicates, mitochondrial reads, and reads that are not properly paired. Reads
are shifted 4 bp to the right and reads on the negative strands should be
shifted 5 bp to the left to accommodate   Tn5 cut positioning. Peaks are then
called with Macs2, and sequence coverage (depth) under peaks is then calculated.

Quality control metrics are also calculated. Trimmed read statistics, mapping
statistics, percent mitochondrial DNA, percent duplication, insert size
estimation, library complexity statistics, and fraction of reads in peak.
