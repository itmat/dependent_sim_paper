input_data_urls = {
    "GSE151923_metaReadCount_ensembl.txt.gz": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151923/suppl/GSE151923%5FmetaReadCount%5Fensembl%2Etxt%2Egz",
    "GSE151923_family.soft.gz": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151923/soft/GSE151923_family.soft.gz",
    "GSE81142.counts.txt.gz": "https://s3.amazonaws.com/itmat.data/tom/dependent_sim_paper/GSE81142.counts.txt.gz",
    "GSE81142_family.soft.gz": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81142/soft/GSE81142_family.soft.gz",
    "BHTC.All_tissues.JTK_only.24h_period.MetaCycle_results.txt.gz": "https://s3.amazonaws.com/itmat.data/BHTC_JTK_RESULTS/BHTC.All_tissues.JTK_only.24h_period.MetaCycle_results.txt.gz",
    "GSE151565_Cortex-counts.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151565/suppl/GSE151565%5FCortex%2Dcounts.csv.gz",
    "plasma_nmr.csv": "https://api2.xialab.ca/api/download/metaboanalyst/plasma_nmr.csv",
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

datasets = {
  "GSE151923": {
        "counts": "data/GSE151923_metaReadCount_ensembl.txt.gz",
        "metadata": "processed/GSE151923_sample_metadata.txt",
        "simulated_controls": expand("simulated_data/Mouse.Cortex.Male.{method}.Control.txt", method = ['pca', 'corpcor', 'wishart', 'indep']),
        "spsimseq": "simulated_data/SPsimSeq.GSE151923.txt",
    },

    "GSE81142": {
        "counts": "data/GSE81142.counts.txt.gz",
        "metadata": "processed/GSE81142_sample_metadata.txt",
        "simulated_controls": expand("simulated_data/Fly.WholeBody.Male.{method}.Control.txt", method = ['pca', 'corpcor', 'wishart', 'indep']),
        "spsimseq": "simulated_data/SPsimSeq.GSE81142.txt",
    },
    "GSE151565": {
        "counts": "processed/Cortex_ZT0-counts.csv",
        "metadata":"data/GSE151565_Cortex-counts.csv.gz", #metadata is in column headers
        "simulated_controls": expand("simulated_data/Cortex_120_normalized_simulated_time_series_{method}.csv", method = ['pca', 'corpcor', 'wishart', 'indep']),
        "spsimseq": "simulated_data/SPsimSeq.GSE151565.txt",
    }
}

wildcard_constraints:
    render_target = "[a-zA-Z0-9]+",
    method = "pca|wishart|indep|corpcor",

rule all:
    input:
        "simulated_data/Mouse.Cortex.Male.pca.Control.txt",
        "simulated_data/Fly.WholeBody.Male.pca.Control.txt",
        "simulated_data/Cortex_120_simulated_time_series_pca.csv",
        "simulated_data/SPsimSeq.GSE81142.txt",
        "simulated_data/SPsimSeq.GSE151923.txt",
        #"time_series_data/Mouse.Cortex.k=2.txt",
        "processed/cyclops/real_data/cyclops_estimated_phaselist.csv",
        expand("processed/cyclops/{method}/batch={batch}/cyclops_estimated_phaselist.csv",
          method = ["indep", "pca", "wishart", "corpcor"], batch=range(0,20)),
        "processed/DE/Mouse.Cortex.Male.pca.fdr.csv",
        "processed/compare_to_real_plot.GSE81142.RDS",
        "processed/compare_to_real_plot.GSE151923.RDS",
        "processed/compare_to_real_plot.GSE151565.RDS",
        "processed/benchmark/results.txt",
        "manuscript/html",
        #"manuscript/docx",
        "manuscript/supplemental/html",

rule generate_sif:
    input:
        "containers/{container}.def"
    output:
        "images/{container}.sif"
    resources:
        mem_mb = 6_000
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
        expand("simulated_data/Mouse.Cortex.Male.{{method}}.{group}{suffix}",
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
        expand("simulated_data/Fly.WholeBody.Male.{{method}}.{group}{suffix}",
            group = ["Case", "Control"],
            suffix = [".txt", ".true_values.txt"])
    resources:
        mem_mb = 6_000
    container:
        "images/dependent_sim.sif",
    script:
        "scripts/simulate_data.fly.R"

rule simulate_fly_spsimseq:
    input:
        "data/GSE81142.counts.txt.gz",
        "processed/GSE81142_sample_metadata.txt",
        sif = ancient("images/spsimseq.sif"),
    output:
        "simulated_data/SPsimSeq.GSE81142.txt",
    resources:
        mem_mb = 36_000
    container:
        "images/spsimseq.sif",
    script:
        "scripts/SPsimSeq.fly.R"

rule simulate_mouse_spsimseq:
    input:
        "processed/GSE151923_sample_metadata.txt",
        "data/GSE151923_metaReadCount_ensembl.txt.gz",
        sif = ancient("images/spsimseq.sif"),
    output:
        "simulated_data/SPsimSeq.GSE151923.txt",
    resources:
        mem_mb = 36_000
    container:
        "images/spsimseq.sif",
    script:
        "scripts/SPsimSeq.mouse.R"

rule simulate_mouse_timeseries_spsimseq:
    input:
        "data/GSE151565_Cortex-counts.csv.gz",
        "processed/Cortex_ZT0-counts.csv",
        sif = ancient("images/spsimseq.sif"),
    output:
        "simulated_data/SPsimSeq.GSE151565.txt",
    resources:
        mem_mb = 60_000
    container:
        "images/spsimseq.sif",
    script:
        "scripts/SPsimSeq.mouse_timeseries.R"

rule sim_de:
    input:
        expand("simulated_data/Mouse.Cortex.Male.{method}.{group}{suffix}",
            method = ["indep", "pca", "wishart", "corpcor"],
            group = ["Case", "Control"],
            suffix = [".txt", ".true_values.txt"]),
        expand("simulated_data/Fly.WholeBody.Male.{method}.{group}{suffix}",
            method = ["indep", "pca", "wishart", "corpcor"],
            group = ["Case", "Control"],
            suffix = [".txt", ".true_values.txt"]),
        sif = "images/dependent_sim.sif",
    output:
        expand("processed/DE/{tissue}.Male.{method}.fdr.csv",
            method = ["indep", "pca", "wishart", "corpcor"],
            tissue = ["Mouse.Cortex", "Fly.WholeBody"]),
    resources:
        mem_mb = 6_000
    container:
        "images/dependent_sim.sif",
    script:
        "scripts/simDE.R"

rule simulate_metabolomics:
    input:
        "data/plasma_nmr.csv",
        sif = "images/dependent_sim.sif",
    output:
        "simulated_data/Plasma_metabolomics.{method}.txt"
    container:
        "images/dependent_sim.sif"
    script:
        "scripts/simulate_metabolomics.R"

rule process_time_series:
    input:
        "data/GSE151565_Cortex-counts.csv.gz",
        ancient("images/dependent_sim.sif"),
    output:
        "processed/mean_per_time_Cortex.csv",
        "processed/Cortex_normalized_time_series_data.csv",
        "processed/Cortex_ZT0-counts.csv",
    container:
        "images/dependent_sim.sif",
    script:
        "scripts/process_time_series.R"

rule simulate_time_series:
    input:
        "processed/Cortex_ZT0-counts.csv",
        "processed/mean_per_time_Cortex.csv",
        "images/dependent_sim.sif",
    output:
        expand("simulated_data/Cortex_120_{datatype}_time_series_{{method}}.csv",
            datatype = ["simulated", "normalized_simulated"])
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
        expression = "simulated_data/Cortex_120_normalized_simulated_time_series_{method}.csv"
    output:
        expression = "processed/cyclops_input/batch={batch}.{method}.csv"
    params:
        batch_size = 4
    script:
        "scripts/make_cyclops_input.py"

rule run_cyclops:
    input:
        expression = "processed/cyclops_input/batch={batch}.{method}.csv",
        seedfile = "processed/cyclic_gene_list.txt",
        sif = "images/cyclops.sif",
    output:
        output = "processed/cyclops/{method}/batch={batch}/cyclops_estimated_phaselist.csv",
    params:
        outdir = "processed/cyclops/{method}/batch={batch}/",
    resources:
        mem_mb = 4000
    container:
        "file://images/cyclops.sif",
    shell:
        "julia /runCYCLOPS.jl --infile {input.expression} --seedfile {input.seedfile} --outdir {params.outdir} --Out_Symbol cyclops --Frac_Var 0.99 --DFrac_Var 0.02"

rule run_cyclops_real_data:
    input:
        expression = "processed/Cortex_normalized_time_series_data.csv",
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
        "julia /runCYCLOPS.jl --infile {input.expression} --seedfile {input.seedfile} --outdir {params.outdir} --Out_Symbol cyclops --Frac_Var 0.99 --DFrac_Var 0.02"

def data_for_dataset(wildcards):
    dataset = datasets[wildcards.dataset]
    return [
        dataset['counts'],
        dataset['metadata'],
        dataset['spsimseq'],
        *dataset['simulated_controls'],
    ]
rule compare_to_real:
    input:
        data_for_dataset,
        sif = "images/quarto.sif",
    output:
        "processed/compare_to_real_plot.{dataset}.RDS"
    container:
        "images/quarto.sif",
    script:
        "scripts/compare_to_real.R"

rule compare_to_real_metabolomics:
    input:
        expand("simulated_data/Plasma_metabolomics.{method}.txt", method = ["indep", "pca", "corpcor", "wishart"]),
        sif = "images/quarto.sif",
    output:
        "processed/compare_to_real_plot_metabolomics.RDS"
    container:
        "images/quarto.sif",
    script:
        "scripts/compare_to_real_metabolomics.R"

rule simulate_for_benchmark:
    input:
        "processed/Cortex_ZT0-counts.csv",
        sif = "images/dependent_sim.sif",
    output:
        "processed/benchmark/{method}/{n_genes}.txt"
    benchmark:
        "processed/benchmark/{method}/{n_genes}.benchmark.txt"
    resources:
        mem_mb = 4_000
    container:
        "images/dependent_sim.sif"
    script:
        "scripts/simulate.benchmark.R"

rule simulate_spsimseq_for_benchmark:
    input:
        "processed/Cortex_ZT0-counts.csv",
        sif = ancient("images/spsimseq.sif"),
    output:
        "processed/benchmark/SPsimSeq/{n_genes}.txt"
    resources:
        mem_mb = 90_000,
    benchmark:
        "processed/benchmark/SPsimSeq/{n_genes}.benchmark.txt"
    container:
        "images/spsimseq.sif"
    script:
        "scripts/SPsimSeq.benchmark.R"

rule simulate_vinecopula_for_benchmark:
    input:
        "processed/Cortex_ZT0-counts.csv",
        sif = ancient("images/spsimseq.sif"),
    output:
        "processed/benchmark/vinecopula/{n_genes}.txt"
    resources:
        mem_mb = 90_000,
    benchmark:
        "processed/benchmark/vinecopula/{n_genes}.benchmark.txt"
    container:
        "images/spsimseq.sif"
    script:
        "scripts/vinecopulib.benchmark.R"

rule simulate_mvrnorm_for_benchmark:
    input:
        "processed/Cortex_ZT0-counts.csv",
        sif = ancient("images/spsimseq.sif"),
    output:
        "processed/benchmark/mvrnorm/{n_genes}.txt"
    resources:
        mem_mb = 60_000,
    benchmark:
        "processed/benchmark/mvrnorm/{n_genes}.benchmark.txt"
    container:
        "images/spsimseq.sif"
    script:
        "scripts/mvrnorm.benchmark.R"

rule benchmark:
    input:
        results = expand("processed/benchmark/{method}/{n_genes}.txt",
            method = ["indep", "pca", "wishart", "corpcor", "SPsimSeq", "vinecopula", "mvrnorm"],
            n_genes = [500, 1000, 2000, 4000, 8000, 16000, 32000],
        ),
        benchmarks = expand("processed/benchmark/{method}/{n_genes}.benchmark.txt",
            method = ["indep", "pca", "wishart", "corpcor", "SPsimSeq", "vinecopula", "mvrnorm"],
            n_genes = [500, 1000, 2000, 4000, 8000, 16000, 32000],
        )
    output:
        "processed/benchmark/results.txt",
    script:
        "scripts/benchmark.R"

rule generate_manuscript:
    input:
        expand("processed/DE/{tissue}.Male.{method}.fdr.csv",
            method = ["indep", "pca", "wishart", "corpcor"],
            tissue = ["Mouse.Cortex", "Fly.WholeBody"]),
        "processed/cyclops/real_data/cyclops_estimated_phaselist.csv",
        expand("processed/cyclops/{method}/batch={batch}/cyclops_estimated_phaselist.csv", 
            method = ["indep", "pca", "wishart", "corpcor"], 
            batch=range(0,20)),
        "processed/compare_to_real_plot.GSE151923.RDS",
        index = "manuscript/paper.qmd",
        sif = "images/quarto.sif",
    output:
        directory("manuscript/{render_target}")
    container:
        "images/quarto.sif"
    shell:
        "quarto render {input.index} --output-dir {wildcards.render_target} --to {wildcards.render_target}"

rule generate_supplemental:
    input:
        "processed/compare_to_real_plot.GSE81142.RDS",
        "processed/compare_to_real_plot.GSE151565.RDS",
        "processed/compare_to_real_plot_metabolomics.RDS",
        "processed/benchmark/results.txt",
        index = "manuscript/supplemental/supplemental.qmd",
        sif = "images/quarto.sif",
    output:
        directory("manuscript/supplemental/{render_target}")
    container:
        "images/quarto.sif"
    shell:
        "quarto render {input.index} --output-dir {wildcards.render_target} --to {wildcards.render_target}"
