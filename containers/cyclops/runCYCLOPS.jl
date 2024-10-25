#! /usr/local/bin/julia
## Adapted from Gang Wu see https://github.com/gangwug/Oslops
## License: GPL >= 2

############################################################
#introduction page for julia language https://docs.julialang.org/en/stable/
#parse the input parameters
using ArgParse     #using the module for command-line argument parsing, see examples at 'https://github.com/carlobaldassi/ArgParse.jl/tree/master/examples'
s = ArgParseSettings()
@add_arg_table s begin
    "--infile", "-i"
        help = "the input expression dataset, csv file (including path)"
        required = true
    "--seedfile", "-s"
        help = "the seed gene list, csv file (including path)"
        required = true
    "--outdir", "-o"
        help = "the directory store ouput csv files"
        required = true
    "--Frac_Var"
    help = "Set Number of Dimensions of SVD to maintain this fraction of variance"
    arg_type = Float64
    default = 0.85
    "--DFrac_Var"
    help = "Set Number of Dimensions of SVD to so that incremetal fraction of variance of var is at least this much"
    arg_type = Float64
    default = 0.03
    "--Seed_MinCV"
    help = "Set the minimal CV for filtering the genes in the seed list, genes with CV below this value will be removed"
    arg_type = Float64
    default = 0.07
    "--Seed_MaxCV"
    help = "Set the maximal CV for filtering the genes in the seed list, genes with CV above this value will be removed"
    arg_type = Float64
    default = 0.75
    "--Seed_Blunt"
    help = "Set the blunt number for removing outliers"
    arg_type = Float64
    default = .99
    "--Seed_MinMean"
    help = "Set the mimimal mean expression cutoff for seed list, genes with mean expression level below this value will be removed"
    arg_type = Float64
    default = 0.2
    "--N_best"
    help = "Number of random initial conditions to try for each optimization"
    arg_type = Int
    default  = 40
    "--Seed_Random"
    help = "the random seed number used for CYCLOPS"
    arg_type = Int
    default  = 12345
    "--Out_Symbol"
    help = "The symbol used in the output file"
end
#the actual argument parsing step, parsed_args is a dictionary type
parsed_args = parse_args(ARGS, s)

infile = string(parsed_args["infile"])
outdir = string(parsed_args["outdir"])
seedfile = string(parsed_args["seedfile"])


############################################################
using StatsBase
using MultivariateStats
using Distributions
addprocs(5)

#change to the directory contains CYCLOPS code and import required modules
include("/CYCLOPS/PNAS_CYCLOPS_PROGRAM_SCRIPTS_REVISION_UPLOAD/CYCLOPS_v6_2a_AutoEncoderModule.jl")
include("/CYCLOPS/PNAS_CYCLOPS_PROGRAM_SCRIPTS_REVISION_UPLOAD/CYCLOPS_v6_2a_CircularStats_U.jl")
include("/CYCLOPS/PNAS_CYCLOPS_PROGRAM_SCRIPTS_REVISION_UPLOAD/CYCLOPS_v6_2a_PreNPostprocessModule.jl")
#include("/CYCLOPS/PNAS_CYCLOPS_PROGRAM_SCRIPTS_REVISION_UPLOAD/CYCLOPS_v6_2a_AutoEncoderModule_multi.jl")
#include("/CYCLOPS/PNAS_CYCLOPS_PROGRAM_SCRIPTS_REVISION_UPLOAD/CYCLOPS_v6_2a_MultiCoreModule_Smooth.jl")
include("/CYCLOPS/PNAS_CYCLOPS_PROGRAM_SCRIPTS_REVISION_UPLOAD/CYCLOPS_v6_2a_Seed.jl")
using CYCLOPS_v6_2a_AutoEncoderModule
using CYCLOPS_v6_2a_PreNPostprocessModule
using CYCLOPS_v6_2a_Seed

############################################################
Seed_MinCV              = parsed_args["Seed_MinCV"]
Seed_MaxCV              = parsed_args["Seed_MaxCV"]
Seed_Blunt              = parsed_args["Seed_Blunt"]
Seed_MinMean            = parsed_args["Seed_MinMean"]
Frac_Var                = parsed_args["Frac_Var"]               # Set Number of Dimensions of SVD to maintain this fraction of variance
DFrac_Var               = parsed_args["DFrac_Var"]              # Set Number of Dimensions of SVD to so that incremetal fraction of variance of var is at least this much
Seed_Random             = parsed_args["Seed_Random"]
N_best                  = parsed_args["N_best"]
Out_Symbol              = parsed_args["Out_Symbol"]

srand(Seed_Random)

###########################################################
#read the seed gene list
fullnonseed_data_BHTC=readcsv(seedfile)
bhtc_seeds=fullnonseed_data_BHTC[2:end ,1]

#get the eigen gene expression profile
fullnonseed_data_merge=readcsv(infile)
fullnonseed_data2=hcat(fullnonseed_data_merge[:,1], fullnonseed_data_merge[:,1], fullnonseed_data_merge)       ## the default get seed function assumes first column=probe, 2nd symbol, 3rd entrez (or just text)
alldata_samples=fullnonseed_data2[1,4:end]
alldata_samples2=reshape(alldata_samples, 1, length(alldata_samples))
seed_symbols_bhtc, seed_data_bhtc                   = getseed(fullnonseed_data2,bhtc_seeds,Seed_MaxCV,Seed_MinCV,Seed_MinMean,Seed_Blunt)
seed_data_bhtc                                  = dispersion!(seed_data_bhtc)
outs_bhtc, norm_seed_data_bhtc = GetEigenGenes(seed_data_bhtc,Frac_Var,DFrac_Var,300)
sample_ids = fullnonseed_data_merge[1,2:end]

###########################################################
# Compute the ordering based off the seed genes
estimated_phaselist, bestnet, global_var_metrics = CYCLOPS_Order(outs_bhtc, norm_seed_data_bhtc, N_best)

###########################################################
prefix = string(outdir, "/", Out_Symbol)
# Output values from Eigen step
#output the dispersed seed gene expression and the eigen expression matrix
seed_out = hcat(seed_symbols_bhtc, seed_data_bhtc)
seed_out = vcat(hcat("geneName", alldata_samples2), seed_out)
writecsv(string(prefix, "_SeedgeneExp.csv"), seed_out)

eigen_out = hcat(
    map(string, fill("eigen_", outs_bhtc), 1:outs_bhtc),
    norm_seed_data_bhtc
)
eigen_out = vcat(hcat("eigenName", alldata_samples2), eigen_out)
writecsv(string(prefix, "_EigengeneExp.csv"), eigen_out)

# Output values from the CYCLOPS_Order step
writecsv(string(prefix, "_estimated_phaselist.csv"), vcat(["ID" "phase"], hcat(sample_ids, estimated_phaselist)))
writecsv(string(prefix, "_global_var_metrics.csv"), hcat(["besterror", "F", "Stat_err"], global_var_metrics))
###########################################################


#write the file parameters out
parafile=string(Out_Symbol, "_para.txt")
f = open(string(outdir, "/", parafile), "a")
for (key, val) in parsed_args 
    write(f, "$key : $val \n")
end
write(f, "\nThe important intermediate parameters in ordering with each eigen cluster is listed below:\n")
close(f)

