Variants where called following the `SARS-CoV-2 workflow`_:
Reads were mapped onto {{ snakemake.config["refgenome_name"] }} with `BBMap`_, and both optical and PCR duplicates were removed with Picard_.
The Freebayes_ haplotype-based variant detector was used to call variants per-sample.
Genotyped variants were filtered with the vcffilter_ by removing any sites with estimated probability of not being polymorphic less than phred 20.
Finally, SnpEff_ was used to predict and report variant effects.
In addition, quality control was performed with FastQC_, Samtools_, and Picard_ and aggregated into an interactive report via MultiQC_.

.. _GATK: https://software.broadinstitute.org/gatk/
.. _BBMap: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. _Picard: https://broadinstitute.github.io/picard
.. _SARS-CoV-2 workflow: https://github.com/avilab/sarscov2
.. _SnpEff: http://snpeff.sourceforge.net
.. _MultiQC: http://multiqc.info/
.. _Samtools: http://samtools.sourceforge.net/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Freebayes: https://github.com/ekg/freebayes
.. _vcffilter: https://github.com/vcflib/vcflib#vcflib
