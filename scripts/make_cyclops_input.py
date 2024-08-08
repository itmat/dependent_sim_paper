import re
import sys

print(sys.version)
import pandas

all_data = pandas.read_csv(snakemake.input.expression, index_col=0)

batch = int(snakemake.wildcards.batch)
batch_size = snakemake.params.batch_size
print(f"Processing batch {repr(batch)} of size {repr(batch_size)}")

sample_nums = [
    int(re.match("ZT[0-9]+_sample([0-9]+)", c).groups()[0]) for c in all_data.columns
]
print(sample_nums)
selected_columns = [
    c
    for c, num in zip(all_data.columns, sample_nums)
    if (num - 1) // batch_size == batch
]

print(selected_columns)

selected_data = all_data[selected_columns]
print(selected_data)

selected_data.to_csv(snakemake.output.expression)
