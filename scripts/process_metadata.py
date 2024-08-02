import re
import gzip

infile = snakemake.input.soft

samples = []
current_sample = None
for line in gzip.open(infile, "rt"):
    m = re.match(r"\^SAMPLE = (GSM[0-9]*)", line)
    if m:
        if current_sample is not None:
            samples.append(current_sample)
        current_sample = {
            "GSM_id": m.groups()[0],
            "description": []
        }

    for attribute, regex in zip(snakemake.params.attribute_names, snakemake.params.attribute_regex):
        m = re.match(regex, line)
        if m:
            current_sample[attribute] = m.groups()[0]

if current_sample not in samples:
    samples.append(current_sample)

## Extract the sample ids
#for sample in samples:
#    match_descr = [x for x in sample['description'] if x.startswith(sample['batch'])]
#    assert len(match_descr) == 1
#    sample['sample_id'] = match_descr[0]

import pandas
df = pandas.DataFrame(samples).drop(columns="description")
df.to_csv(snakemake.output.table, sep="\t", index=False)
