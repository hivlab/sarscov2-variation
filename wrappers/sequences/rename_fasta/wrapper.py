from Bio import SeqIO
from hashlib import blake2s

sample = snakemake.params.get("sample", "")
stub = snakemake.params.get("stub", "")
assert len(sample) > 0, "Sample name is missing"
h = blake2s()

with open(snakemake.input[0], "r") as input_handle, open(
    snakemake.output[0], "w"
) as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        h.update(sample.encode("utf-8"))
        record.id = stub.format(h.hexdigest()[:10])
        record.description = ""
        SeqIO.write(record, output_handle, "fasta")
