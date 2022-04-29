# Allele-specific ATAC-seq workflow

Nextflow workflow to process raw ATAC-seq data for allele-specific analysis. It takes compressed fastq files as input and generates quality metrics, peak files based on total signal ("Gall") using 3 different approaches (macs2: fragments, cut sites; HMMRATAC) and normalized signal tracks of fragments, cut sites and nucleosome-free regions. Also, it saves log files that can be parsed to calculate the number of reads filtered out in each step.

## Author

Yuvia Alhelí PÉREZ-RICO

Affiliation: European Molecular Biology Laboratory | EMBL

## Dependencies

To use this workflow, you will need to install the following programs, the indicated version is the one that I have used to analyze my data:

- nextflow (20.04.1)
- FastQC (v0.11.9)
- Trim Galore (0.6.3)
- Bowtie2 (2.3.4.1)
- samtools (1.9)
- SNPsplit (0.3.4)
- Picard Tools (2.20.8)
- deepTools (3.4.3)
- bedtools (v2.29.2)
- macs2 (2.2.7.1)
- HMMRATAC (1.2.10)
- bedGraphToBigWig (v 4)
- R (4.0.2)
- ATACseqQC bioconductor package (1.12.0)
- bioconductor packages for peak annotation: ggupset (0.3.0), org.Mm.eg.db (3.11.4), ChIPseeker (1.24.0) and TxDb.Mmusculus.UCSC.mm10.knownGene (3.10.0).
- ENCODE libraries: encode_lib_common.py and encode_lib_genomic.py, part of the ENCODE-DCC atac-seq-pipeline
- ENCODE script: encode_task_tss_enrich.py, part of the ENCODE-DCC atac-seq-pipeline

Note: I suggest to install all programs in a conda environment, except the packages used to annotate peaks and test for signal enrichment at TSSs because I experienced dependency issues. For instance, python2 is required to run the ENCODE script. Thus, I use independent environments for those processes (lines 797 and 841 of the main script).

## How to use the workflow

Thank you for your interest in using the workflow!

After downloading the main script, the configuration file and additional scripts within the bin folder, you will need to prepare the following files and change the paths in the configuration file accordingly:

- N-masked bowtie2 index using SNPs of 2 strains.
- Compressed file with the SNPs of the 2 strains for SNPsplit.
- Tab-separated file with chromosome sizes (only 2 columns: chromosome name and size).
- Bed file containing blacklist genomic regions.
- Bed file containing TSS annotations.
- Comma-separated file indicating the sample name and paths to the fastq files. An example is available in the 'docs' folder.

If you are not using SLURM, then do additional modifications to the configuration file considering the workload manager that you use. Finally, write a simple bash script that will be submitted to the cluster to activate the conda environment and start the main nextflow job:

`source /home/user/miniconda2/bin/activate /home/user/conda-envs/ATAC-seq`

`nextflow run allelic_ATAC-seq.nf`

Please, keep in mind that this workflow was written for the processing of mouse paired-end data, therefore, you will need to do some changes to the main and additional scripts before using it with data of other species, for example, changing the strain names and databases for peak annotation. Also, I skip some chromosomes for the normalization of signal tracks, as I know that they tend to get duplicated in the cell lines that are commonly handled in our lab or to ease the visualization of differentiated and undifferentiated cells. Modify the 'signal_tracks1' and 'signal_tracks2' processes in case that you do not want to skip chromosomes.

## Acknowledgements

This repository is part of a project that has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 882771.

## License

Licenced under the GNU general public license.

