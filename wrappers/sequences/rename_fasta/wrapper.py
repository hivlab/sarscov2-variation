from snakemake.shell import shell
from hashlib import blake2s

sample = snakemake.params.get("sample", "")
stub = snakemake.params.get("stub", "")
assert len(sample) > 0, "Sample name is missing"
h = blake2s()
h.update(sample.encode("utf-8"))
new_id = stub + h.hexdigest()

shell(
    "sed 's/^>.*$/>{new_id}/g' {snakemake.input[0]} > {snakemake.output[0]}"
    )
