from Bio import Entrez
import os

Entrez.email = snakemake.params["email"]
Entrez.api_key = snakemake.params["api_key"]
handle = Entrez.esearch(db="nucleotide", retmax=5000, term="SARS-CoV-2", idtype="acc")
record = Entrez.read(handle)
handle.close()
acc = record["IdList"]
handle = Entrez.efetch(db="nucleotide", id=",".join(acc), rettype="gb", retmode="text")
with open(snakemake.output[0], "w") as h:
    h.writelines(handle)
handle.close()
