import pathlib
import re
import pandas

dataset = {}
for filepath in pathlib.Path("data/GSE81142_RAW").glob("*.txt.gz"):
    m = re.match(r"GSM[0-9]*_([a-zA-Z0-9]*).reverse.htseq.counts.txt.gz", filepath.name)
    if m:
        sample_name = m.groups()[0]
    else:
        print(f"WARNING: Skipping unexpected file name {filepath}")
        continue

    d = pandas.read_csv(filepath, sep="\t", names=["GENE_ID", "counts"])
    d = d[~d.GENE_ID.str.startswith("__")]
    dataset[sample_name] = d.set_index("GENE_ID").counts

pandas.DataFrame(dataset).to_csv("processed/GSE81142.counts.txt.gz", sep="\t")

metadata = []
for sample in dataset.keys():
    m = re.match(f"M([mf])([0-9])h([0-9])cR([0-9])", sample)
    if not m:
        print(f"ERROR: sample name not of the expected format: {sample}")
        break
    sex,hour,concentration,rep = m.groups()
    metadata.append({
        "sample_id": sample,
        "sex": {'m': 'male', 'f': 'female'}[sex],
        "hour": int(hour),
        "concentration": int(concentration),
        "rep": int(rep),
    })
df = pandas.DataFrame(metadata)
df.to_csv("processed/GSE81142.metadata.txt", sep="\t", index=False)
