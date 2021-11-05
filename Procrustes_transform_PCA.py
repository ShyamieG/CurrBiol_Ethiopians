## SCRIPT: Procrustes_transform_PCA.py
## AUTHOR: SSG
## OVERVIEW: Produces a set of PCAs by dropping in a set of individuals one at a time
## REQUIRED ARGUMENTS:
##         --bfile (stem of the PLINK binary file name that contains data on both reference individuals and those to be dropped)
##         --sample_list (a list of samples to drop in, having at least the first two columns of the .fam file)
## OPTIONAL ARGUMENTS:
##         --drop_missing_data (should SNPs not present in the sample being dropped in be removed when doing the PCA? default True)
##         --make_pseudohap (should reference samples be made pseudohaploid, default False)
##         --n_eigs (number of eigenvalues to compute, default 6)
##         --LD_winsize (window size in SNPs - to pass to PLINK for indep-pairwise LD filter)
##         --LD_stepsize (step size in SNPs - to pass to PLINK for indep-pairwise LD filter)
##         --LD_r2 (maximum pairwise SNP correlation to tolerate  - to pass for indep-pairwise LD filter)
##         --maf_min (minimum minor allele frequency - to pass to PLINK for maf filter)
##
## USAGE EXAMPLES:
## $ python Procrustes_transform_PCA.py --bfile PCA_data --sample_list drop_in.list --n_eigs 10 

#! /usr/bin/env python
import argparse as ap
import os
import sys

import pandas as pd
import numpy as np

### Define main function
def main(args):
    ## Parse arguments
    argp = ap.ArgumentParser(description="")
    argp.add_argument("--bfile", required=True)
    argp.add_argument("--sample_list", required=True)
    argp.add_argument("--drop_missing_data", default=True, type=bool)
    argp.add_argument("--make_pseudohap", default=False, type=bool)
    argp.add_argument("--n_eigs", default=6, type=int)
    argp.add_argument("--LD_winsize", type=int, required=False)
    argp.add_argument("--LD_stepsize", type=int, required=False)
    argp.add_argument("--LD_r2", type=float, required=False)
    argp.add_argument("--maf_min", type=float, required=False)
    args = argp.parse_args(args)

    ## Perform a PCA for each dropped in individual
    print("Parsing input...")
    all_samples = pd.read_csv(str(args.bfile) + ".fam", sep="\s+", header=None)
    dropin_samples = pd.read_csv(str(args.sample_list), sep="\s+", header=None)
    idx_ref = []
    dropin_samp_concatIDs = (dropin_samples.loc[0:,0] + dropin_samples.loc[0:,1]).values.tolist()
    all_samp_concatIDs = (all_samples.loc[0:,0] + all_samples.loc[0:,1]).values.tolist()
    for i in all_samp_concatIDs:
        if i not in dropin_samp_concatIDs:
            idx_ref.append(all_samp_concatIDs.index(i))

    ref_samples = all_samples.iloc[idx_ref,]
    
    if args.make_pseudohap == True:
        print("Making all samples pseudo-haploid...")
        import pandas_plink as pdp
        import random
        bedfile = pdp.read_plink1_bin(str(args.bfile) + ".bed")
        nInds = bedfile.shape[0]
        for ind_idx in range(nInds):
            print("    working on sample " + str(bedfile.iid.values[ind_idx]) + " (" + str(ind_idx+1) + "/" + str(nInds) + ")")
            genos = bedfile.sel(sample=bedfile.sample.values[ind_idx]).values.tolist()
            het_idx = np.where(np.array(genos)==1)[0][0:]
            for g in het_idx:
                genos[g] = float(np.array(random.sample([0., 2.], k=1)))
            alleles = [0] * len(genos)
            for n in range(len(genos)):
                if genos[n]==0.:
                    alleles[n] = bedfile.a0.values[n]
                elif genos[n]==2.:
                    alleles[n] = bedfile.a1.values[n]
            IND_ped = [bedfile.fid.values[ind_idx].tolist()] + [bedfile.iid.values[ind_idx].tolist()] + [bedfile.father.values[ind_idx].tolist()] + [bedfile.mother.values[ind_idx].tolist()] + [bedfile.gender.values[ind_idx].tolist()] + [str(bedfile.trait.values[ind_idx].tolist())] + [x for pair in zip(alleles, alleles) for x in pair]
            if ind_idx==0:
                with open(str(args.bfile) + "_pseudoHap.ped", "w") as pedfile:
                    pedfile.write("%s" % " ".join(map(str, IND_ped)))
            else:
                with open(str(args.bfile) + "_pseudoHap.ped", "a") as pedfile:
                    pedfile.write("\n" + "%s" % " ".join(map(str, IND_ped)))
        os.system("awk \'{print $1, $2, $3, $4}\' " + str(args.bfile)+ ".bim > " + str(args.bfile) + "_pseudoHap.map")
        args.bfile = str(args.bfile) + "_pseudoHap"
        os.system("plink2 --file " + str(args.bfile) + " --make-bed --out " + str(args.bfile))
        os.system("rm " + str(args.bfile) + ".ped")
        os.system("rm " + str(args.bfile)+ ".map")

    # Modify to run in parallel?
    os.system("echo sample_ID n_SNPs > dropInSamplesInfo")
    print("Doing individual PCAs...")
    for i in range(len(dropin_samples)):
        # Make files for individual i
        IND_df = dropin_samples.iloc[i:i+1,]
        ind = IND_df.iloc[0,1]
        print("    working on sample " + str(ind) + "...")
        rem_df = dropin_samples.drop(i)
        rem_df.to_csv(".REM.list", sep=" ", index=False, header=False)
        if args.drop_missing_data == True:
            # Keep individual i, and only the SNPs typed in that individual
            IND_df.to_csv(".KEEP.list", sep=" ", index=False, header=False)
            os.system("plink2 --bfile " + str(args.bfile) + " --keep .KEEP.list --geno 0 --make-bed --out .KEEP")
            # Remove everyone else that needs to be dropped in, and keep only SNPs typed in individual i
            os.system("plink2 --bfile " + str(args.bfile) + " --remove .REM.list --extract .KEEP.bim --make-bed --out .TEMP")
        else:
            os.system("plink2 --bfile " + str(args.bfile) + " --remove .REM.list --make-bed --out .TEMP")
        
        if (args.LD_winsize is not None) & (args.LD_stepsize is not None) & (args.LD_r2 is not None):
            os.system("plink2 --bfile .TEMP --indep-pairwise " + str(args.LD_winsize) + " " + str(args.LD_stepsize) + " " + str(args.LD_r2) + " --out .LD_filter")
            os.system("plink2 --bfile .TEMP --extract .LD_filter.prune.in --make-bed --out .TEMP2")
            os.system("mv .TEMP2.fam .TEMP.fam"); os.system("mv .TEMP2.bim .TEMP.bim"); os.system("mv .TEMP2.bed .TEMP.bed"); os.system("rm .TEMP2.*")
        
        if args.maf_min is not None:
            os.system("plink2 --bfile .TEMP --maf " + str(args.maf_min) + " --make-bed --out .TEMP2")
            os.system("mv .TEMP2.fam .TEMP.fam"); os.system("mv .TEMP2.bim .TEMP.bim"); os.system("mv .TEMP2.bed .TEMP.bed"); os.system("rm .TEMP2.*")
         
        # Create parameter file for smartpca, and run
        os.system("echo genotypename: .TEMP.bed > .par.all")
        os.system("echo snpname: .TEMP.bim >> .par.all")
        os.system("echo indivname: .TEMP.fam >> .par.all")
        os.system("echo output: " + str(ind) + "_PCA.eigs >> .par.all")
        os.system("echo altnormstyle: NO >> .par.all")
        os.system("echo numeigs: " + str(args.n_eigs) + " >> .par.all")
        os.system("echo numoutlieriter: 0 >> .par.all")
        os.system("echo numoutlierevec: 0 >> .par.all")
        os.system("echo missingmode: NO >> .par.all")
        os.system("echo nsnpldregress: 0 >> .par.all")
        os.system("echo noxdata: YES >> .par.all")
        os.system("echo nomalexhet: YES >> .par.all")
        os.system("smartpca -p .par.all")
        
        # Record SNP number for each drop in sample
        os.system("echo " + str(ind) + " `wc -l .TEMP.bim | awk '{print $1}'` >> dropInSamplesInfo")
    
                
    ## Clean up
    os.system("rm .KEEP.* .REM.list .TEMP.* .par.all .LD_filter.*")

    ## Run Procrustes transformation with R script, package "vegan"
    print("Running Procrustes transformation...")
    os.system("Rscript Procrustes_transform_PCA.R " + str(args.bfile) + " " + str(args.n_eigs))

    return 0

### Run script
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

