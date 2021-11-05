## SCRIPT: SNP_array_match.py
## AUTHOR: SSG
## OVERVIEW: Downsample data to match a given SNP array in SNP density and minor allele frequency distribution
## This script accepts two binary PLINK file names - one for the data to be downsampled, and for the template SNP array
## It returns a list of SNPs to extract (i.e. with PLINK) from the data file
## This script goes along each 1Mb window of the reference SNP array and selects a similar* number of SNPs in the target data that also roughly matches the allele frequency distribution of the reference SNP array for that window

## REQUIRED ARGUMENTS:
##         --data (the data to be downsampled, binary PLINK format)
##         --target (target SNP array data to match against, binary PLINK format)
##         --outdir directory to output the file of variants to 
## OPTIONAL ARGUMENTS:
##         --remove_mono (remove all monomorphic SNPs from target)
##         --norm_target (adjust base pair positions of target variants so each chromosome starts at 0)
##
## USAGE EXAMPLE:
## python SNP_array_match.py --data steadyNe_n100 --target MAJANG_template --outdir /share/hennlab/projects/Ethiopians/ --drop_target_mono True --norm_target True --add_mono_vars False --fold_target True

#! /usr/bin/env python
import argparse as ap
import os
import sys

import pandas as pd
import numpy as np
import random

### Define main function
def main(args):
    ## Parse arguments
    argp = ap.ArgumentParser(description="")
    argp.add_argument("--data", required=True)
    argp.add_argument("--target", required=True)
    argp.add_argument("--drop_target_mono", default=True, type=bool)
    argp.add_argument("--norm_target", default=True, type=bool)
    argp.add_argument("--add_mono_vars", default=True, type=bool)
    argp.add_argument("--fold_target", default=True, type=bool)
    argp.add_argument("--outdir", required=True)
    args = argp.parse_args(args)

    ## Read .bim and .afreq files in
    print("Reading in data")
    DATA_bim = pd.read_csv(args.data + ".bim", header=None, sep="\t", names=["chr", "varID", "gen_pos", "bp_pos", "A1", "A2"])
    DATA_freq = pd.read_csv(args.data + "_freq.frq", sep="\t", header=0, names=["chr", "varID", "A1", "A2", "freq", "count"] )
    DATA = pd.concat([DATA_bim[["chr", "varID", "bp_pos"]], DATA_freq[["freq"]]], axis=1)
    TARGET_bim = pd.read_csv(args.target + ".bim", header=None, sep ="\t", names=["chr", "varID", "gen_pos", "bp_pos", "A1", "A2"])
    TARGET_freq = pd.read_csv(args.target + "_freq.frq", sep="\t", header=0, names=["chr", "varID", "A1", "A2", "freq", "count"])
    if args.fold_target==True:
        TARGET_freq.loc[TARGET_freq["freq"] > 0.5, "freq"] = 1-TARGET_freq.loc[TARGET_freq["freq"] > 0.5, "freq"]

    TARGET = pd.concat([TARGET_bim[["chr", "varID", "bp_pos"]], TARGET_freq[["freq"]]], axis=1)
    print(TARGET)
    # Clean up a bit
    #os.system("rm -f " + args.data + "_freq.*")
    #os.system("rm -f " + args.target + "_freq.*")

    ## 'Normalize' target base pair positions, if specified
    if args.norm_target==True:
        print("Normalizing target array")
        for CHR in np.unique(np.array(TARGET["chr"])):
            TARGET.loc[(TARGET.chr==CHR), "bp_pos"] = TARGET[TARGET.chr==CHR]["bp_pos"] - min(TARGET[TARGET.chr==CHR]["bp_pos"])

    ## Exclude monomorphic SNPs from target set
    if args.drop_target_mono == True:
        mono_vars = TARGET[(TARGET.freq == 0) | (TARGET.freq == 1)]
        print("Removing monomorphic variants from target")
        TARGET.drop(labels=list(mono_vars.index), axis=0, inplace=True)

    ## Sample variants in the data that match the target
    selected_vars = pd.DataFrame()
    for CHR in np.unique(np.array(TARGET["chr"])):
        print("Sampling chromosome " + str(CHR), end="... "); sys.stdout.flush()
        # Calculate SNP density
        Density = np.histogram(TARGET[TARGET.chr==CHR]["bp_pos"], bins=range(0, max(TARGET[TARGET.chr==CHR]["bp_pos"])+1000000, 1000000))
        Data_Dens = np.histogram(DATA[DATA.chr==CHR]["bp_pos"], bins=range(0, max(TARGET[TARGET.chr==CHR]["bp_pos"])+1000000, 1000000))
        for i in list(Density[1][1:]):
            idx=list(Density[1][1:]).index(i)
            if Density[0][idx]>0:
                # Calculate allele frequency distribution
                freqs = TARGET.loc[(TARGET.chr==CHR) & (i-1000000 < TARGET.bp_pos) & (TARGET.bp_pos < i), "freq"]
                freq_density = np.histogram(freqs, bins=np.arange(0, 1, 0.05))
                DATA_region = DATA.loc[(DATA.chr==CHR) & (DATA.bp_pos > i-1000000) & (DATA.bp_pos < i)] # this is the problem
                for j in list(freq_density[1][1:]):
                    jdx = list(freq_density[1][1:]).index(j)
                    if freq_density[0][jdx] > 0:
                        # Identify all variants in DATA that match
                        matching_vars = list(pd.DataFrame(DATA_region.loc[(DATA_region.freq > j-0.05) & (DATA_region.freq <= j), "varID"])["varID"])
                       # print(matching_vars)
                        # If more elements to choose from than need, randomly sample them
                        if len(matching_vars) > freq_density[0][jdx]:
                            sampled_vars = random.sample(matching_vars, freq_density[0][jdx])
                            selected_vars = selected_vars.append(sampled_vars, ignore_index=True)
                        # Otherwise, take all of them
                        elif len(matching_vars) > 0 & len(matching_vars) <= freq_density[0][jdx]:
                            selected_vars = selected_vars.append(matching_vars, ignore_index=True)
        print("done")

    ## Put the selected variants into the original order
    
    selected_vars_sorted = (DATA.varID[DATA.varID.isin(selected_vars[0])].unique())
    
    OUT = args.outdir + "/" + args.data + "_vars_to_extract.list" 
    ## Write out a list of variants to a file
    selected_vars_sorted.tofile(OUT, sep="\n", format="%s")

    ## Extract these from the original dataset
#    print("Extracting selected variants from original data")
#    os.system("plink2 --bfile " + args.data + " --extract vars_to_extract.list --make-bed --out " + args.data + "_SNP-arrayThinned --silent")
#    os.system("rm -f vars_to_extract.list")

    ## Add monomorphic SNPs back in
#    if args.add_mono_vars==True:`
#        print("Adding monomorphic variants to data")
#        mono_vars.index = np.arange(len(mono_vars))
#        fam = pd.read_table(args.data + ".fam", header=None)
#        ped = pd.concat([fam, pd.DataFrame(np.full(shape=(len(fam), len(mono_vars)*2), fill_value="."))], axis=1)
#        map = pd.concat([mono_vars.chr, mono_vars["chr"].astype(str) + ":" + mono_vars["bp_pos"].astype(str), pd.Series("0", range(len(mono_vars))), mono_vars.bp_pos], axis=1)
#        print("Writing monomorphic variant PLINK files")
#        ped.to_csv("mono_vars.ped", sep="\t", index=False, header=False)
#        map.to_csv("mono_vars.map", sep="\t", index=False, header=False)
#        print("Converting to binary files")
#        os.system("plink2 --file mono_vars --make-bed --out mono_vars --silent")
#        os.system("plink2 --bfile mono_vars --bmerge " + args.data + "_SNP-arrayThinned --make-bed --out " + args.data + "_SNP-arrayThinned_monoVarsAdded --silent")
#        os.system("rm -f mono_vars.*")

    return 0

### Run script
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

