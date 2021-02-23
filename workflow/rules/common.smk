
def get_fastq(wildcards):
    data = pep.sample_table
    read_cols = [col for col in data if "path" in col]
    data.set_index(["run"], append=True, inplace=True, drop=False)
    reads = data.loc[(wildcards.sample, wildcards.run), read_cols].dropna()
    assert len(read_cols) in [1, 2], "Enter one or two FASTQ file paths"
    if len(read_cols) == 2:
        return {"in1": reads["read1_path"], "in2": reads["read2_path"]}
    else:
        return {"input": reads["read1_path"]}
