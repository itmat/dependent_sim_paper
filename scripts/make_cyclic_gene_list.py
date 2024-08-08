import polars as pl

bhtc = pl.read_csv("data/BHTC.All_tissues.JTK_only.24h_period.MetaCycle_results.txt.gz", separator="\t", null_values=["NaN", "NA"])

frac_tissues_sig = bhtc.group_by(["ensID", "geneSymbol"]) \
                       .agg((pl.col("JTK_pvalue") < 0.05).mean().alias("frac_sig"))

# Take genes that were rhythmic in at least 75% of the tissues
valid = frac_tissues_sig \
        .filter(pl.col("frac_sig") >= 0.75) \
        .select("ensID", "geneSymbol")

# Some entries have multiple IDs separated by a ','
# we'll use all
valid_expanded = valid.select(
        pl.col("ensID")
            .str.split(",")
            .explode()
    )\
    .filter(pl.col("ensID").is_not_null())

valid_expanded.write_csv("processed/cyclic_gene_list.txt", separator="\t")
