from Bio import SeqIO
from hashlib import blake2s

h = blake2s()

with open("output/consensus_fix.fa", "r") as input_handle, open("output/consensus_ano.fa", "w") as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        h.update(record.id.encode("utf-8"))
        record.id = h.hexdigest()
        record.description = h.hexdigest()
        SeqIO.write(record, output_handle, "fasta")