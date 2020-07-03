Variants where called following the `SARS-CoV-2 workflow`_:
Reads were mapped onto {{ snakemake.config["refgenome_name"] }} with `BBMap`_, and both optical and PCR duplicates were removed with `Clumpify`_, reads were quality trimmed with `BBDuk`_.
The Lofreq_ variant detector was used to call variants per-sample.
Genotyped variants were filtered with the vcffilter_ by removing any sites with estimated probability of not being polymorphic less than phred 30 and allele frequency less than 0.5.
Finally, SnpEff_ was used to predict and report variant effects.
In addition, quality control was performed with FastQC_, Samtools_, and fastq_screen_ and aggregated into an interactive report via MultiQC_.
Consensus sequence was generated using Samtools_, nucleotide poisitions with coverage less than 20 were masked.

.. _BBMap: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. _Clumpify: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. _BBDuk: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. _SARS-CoV-2 workflow: https://github.com/avilab/sarscov2
.. _SnpEff: http://snpeff.sourceforge.net
.. _MultiQC: http://multiqc.info/
.. _Samtools: http://samtools.sourceforge.net/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Lofreq: https://csb5.github.io/lofreq/
.. _vcffilter: https://github.com/vcflib/vcflib#vcflib
.. _fastq_screen: https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
