input_data_urls = {
    "GSE151923_metaReadCount_ensembl.txt.gz": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151923/suppl/GSE151923%5FmetaReadCount%5Fensembl%2Etxt%2Egz",
    "GSE151923_family.soft.gz": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151923/soft/GSE151923_family.soft.gz",
    "GSE81142.counts.txt.gz": "https://s3.amazonaws.com/itmat.data/tom/dependent_sim_paper/GSE81142.counts.txt.gz",
    "GSE81142_family.soft.gz": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81142/soft/GSE81142_family.soft.gz",
    "BHTC.All_tissues.JTK_only.24h_period.MetaCycle_results.txt.gz": "https://s3.amazonaws.com/itmat.data/BHTC_JTK_RESULTS/BHTC.All_tissues.JTK_only.24h_period.MetaCycle_results.txt.gz",
    "GSE151565_Liver-counts.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151565/suppl/GSE151565%5FLiver%2Dcounts%2Ecsv%2Egz",
}
soft_metadata_attributes = {
  "GSE151923": {
    "title": r"!Sample_title = (.*)",
    "geo_accesion": r"!Sample_geo_accession = (.*)",
    "sample_id": r"!Sample_description = ([a-zA-Z0-9]+_[a-zA-Z0-9]+_[a-zA-Z0-9]+)",
    "sex": r"!Sample_characteristics_ch1 = Sex: (.*)",
    "batch": r"!Sample_characteristics_ch1 = batch: (.*) group",
  },
  "GSE81142": {
    "title": r"!Sample_title = (.*)",
    "geo_accesion": r"!Sample_geo_accession = (.*)",
    "sex": r"!Sample_characteristics_ch1 = Sex: (.*)",
    "concentration": r"!Sample_characteristics_ch1 = dmso % v/v: ([\.0-9]+)",
    "hour": r"!Sample_characteristics_ch1 = time (h): ([0-9]+)",
  },
}

rule all:
    input:
        "simulated_data/Mouse.Cortex.Male.k=2.Control.txt",
        "simulated_data/Fly.WholeBody.Male.k=2.Control.txt",
        "simulated_data/Liver_120_simulated_time_series_k=2.csv",
        #"time_series_data/Mouse.Cortex.k=2.txt",
        "processed/cyclops/real_data/cyclops_estimated_phaselist.csv",
        expand("processed/cyclops/k={k}/batch={batch}/cyclops_estimated_phaselist.csv",
          k=[0,2], batch=range(0,20)),

rule generate_sif:
    input:
        "containers/{container}.def"
    output:
        "images/{container}.sif"
    shell:
        "apptainer build {output} {input}"

rule download_data:
  output:
    "data/{data_file}"
  params:
    url = lambda wildcards: input_data_urls[wildcards.data_file]
  shell:
    "curl -o {output} {params.url} "

rule process_soft_metadata:
  input:
    soft = "data/{GSE}_family.soft.gz"
  output:
    table = "processed/{GSE}_sample_metadata.txt"
  params:
    attribute_names = lambda wildcards: list(soft_metadata_attributes[wildcards.GSE].keys()),
    attribute_regex = lambda wildcards: list(soft_metadata_attributes[wildcards.GSE].values()),
  script:
    "scripts/process_metadata.py"

rule simulate_mouse:
    input:
        "processed/GSE151923_sample_metadata.txt",
        "data/GSE151923_metaReadCount_ensembl.txt.gz",
        sif = "images/dependent_sim.sif",
    output:
        expand("simulated_data/Mouse.Cortex.Male.k={k}.{group}{suffix}",
            k = [0,2],
            group = ["Case", "Control"],
            suffix = [".txt", ".true_values.txt"])
    resources:
        mem_mb = 6_000
    container:
        "images/dependent_sim.sif",
    script:
        "scripts/simulate_data.R"

rule simulate_fly:
    input:
        "data/GSE81142.counts.txt.gz",
        "processed/GSE81142_sample_metadata.txt",
        sif = "images/dependent_sim.sif",
    output:
        expand("simulated_data/Fly.WholeBody.Male.k={k}.{group}{suffix}",
            k = [0,2],
            group = ["Case", "Control"],
            suffix = [".txt", ".true_values.txt"])
    resources:
        mem_mb = 6_000
    container:
        "images/dependent_sim.sif",
    script:
        "scripts/simulate_data.fly.R"

rule process_time_series:
    input:
        "data/GSE151565_Liver-counts.csv.gz",
        "images/dependent_sim.sif",
    output:
        "processed/mean_per_time_Liver.csv",
        "processed/Liver_normalized_time_series_data.csv",
        "processed/Liver_ZT0-counts.csv",
    container:
        "images/dependent_sim.sif",
    script:
        "scripts/process_time_series.R"

rule simulate_time_series:
    input:
        "processed/Liver_ZT0-counts.csv",
        "processed/mean_per_time_Liver.csv",
        "processed/GSE151923_sample_metadata.txt",
        "data/GSE151923_metaReadCount_ensembl.txt.gz",
        "images/dependent_sim.sif",
    output:
        "simulated_data/Liver_120_simulated_time_series_k=2.csv",
        "simulated_data/Liver_120_normalized_simulated_time_series_k=2.csv",
        "simulated_data/Liver_120_simulated_time_series_k=0.csv",
        "simulated_data/Liver_120_normalized_simulated_time_series_k=0.csv",
    resources:
        mem_mb = 6_000
    container:
        "images/dependent_sim.sif"
    script:
        "scripts/simulate_time_data.R"

rule make_cyclic_gene_list:
    input:
        "data/BHTC.All_tissues.JTK_only.24h_period.MetaCycle_results.txt.gz"
    output:
        "processed/cyclic_gene_list.txt"
    script:
        "scripts/make_cyclic_gene_list.py"

rule make_cyclops_input:
    input:
        expression = "simulated_data/Liver_120_normalized_simulated_time_series_k={k}.csv"
    output:
        expression = "processed/cyclops_input/batch={batch}.k={k}.csv"
    params:
        batch_size = 6
    script:
        "scripts/make_cyclops_input.py"

rule run_cyclops:
    input:
        expression = "processed/cyclops_input/batch={batch}.k={k}.csv",
        seedfile = "processed/cyclic_gene_list.txt",
        sif = "images/cyclops.sif",
    output:
        output = "processed/cyclops/k={k}/batch={batch}/cyclops_estimated_phaselist.csv",
    params:
        outdir = "processed/cyclops/k={k}/batch={batch}/",
    container:
        "file://images/cyclops.sif",
    shell:
        "julia /runCYCLOPS.jl --infile {input.expression} --seedfile {input.seedfile} --outdir {params.outdir} --Out_Symbol cyclops"

rule run_cyclops_real_data:
    input:
        expression = "processed/Liver_normalized_time_series_data.csv",
        seedfile = "processed/cyclic_gene_list.txt",
        sif = "images/cyclops.sif",
    output:
        output = "processed/cyclops/real_data/cyclops_estimated_phaselist.csv",
    params:
        outdir = "processed/cyclops/real_data/",
    container:
        "file://images/cyclops.sif",
    resources:
        mem_mb = 6_000,
    shell:
        "julia /runCYCLOPS.jl --infile {input.expression} --seedfile {input.seedfile} --outdir {params.outdir} --Out_Symbol cyclops"
