Variants where called following the `SARS-CoV-2 workflow`_:
Reads were mapped onto {{ snakemake.config["refgenome_name"] }} with `BBMap`_, and both optical and PCR duplicates were removed with Picard_, followed by base recalibration with GATK_.
The GATK_ HaplotypeCaller was used to call variants per-sample, including summarized evidence for non-variant sites (GVCF_ approach).
Then, GATK_ genotyping was done in a joint way over GVCF_ files of all samples.
{% if snakemake.config["filtering"]["vqsr"] %}
Genotyped variants were filtered with the GATK_ VariantRecalibrator approach.
{% else %}
Genotyped variants were filtered using hard thresholds.
For SNVs, the criterion ``{{ snakemake.config["filtering"]["hard"]["snvs"] }}`` was used, for Indels the criterion ``{{ snakemake.config["filtering"]["hard"]["indels"] }}`` was used.
{% endif %}
Finally, SnpEff_ was used to predict and report variant effects.
In addition, quality control was performed with FastQC_, Samtools_, and Picard_ and aggregated into an interactive report via MultiQC_.

.. _GATK: https://software.broadinstitute.org/gatk/
.. _BBMapP: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. _Picard: https://broadinstitute.github.io/picard
.. _SARS-CoV-2 workflow: https://github.com/avilab/sarscov2
.. _GVCF: https://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf
.. _SnpEff: http://snpeff.sourceforge.net
.. _MultiQC: http://multiqc.info/
.. _Samtools: http://samtools.sourceforge.net/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
