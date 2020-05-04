from Bio import SeqIO
from hashlib import blake2s

sample = snakemake.params.get("sample", "")
stub = snakemake.params.get("stub", "")
assert len(sample) > 0, "Sample name is missing"
h = blake2s()

with open(snakemake.input, "r") as input_handle, open(
    snakemake.output, "w"
) as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        h.update(sample.encode("utf-8"))
        record.id = h.hexdigest()
        record.description = stub + h.hexdigest()
        SeqIO.write(record, output_handle, "fasta")
