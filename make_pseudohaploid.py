## SCRIPT: make_pseudohaploid.py
## AUTHOR: SSG
## OVERVIEW: Forces sample to be pseudohaploid by randomly picking one allele at heterozygous positions
## REQUIRED ARGUMENTS:
##         --bfile (stem of the PLINK binary file name that contains data on both reference individuals and those to be dropped)
##
## USAGE EXAMPLES:
## $ python make_pseudohap.py --bfile PCA_data

#! /usr/bin/env python
import argparse as ap
import os
import sys

import pandas_plink as pdp
import random
import numpy as np

### Define main function
def main(args):
    ## Parse arguments
    argp = ap.ArgumentParser(description="")
    argp.add_argument("--bfile", required=True)
    args = argp.parse_args(args)

    bedfile = pdp.read_plink1_bin(str(args.bfile) + ".bed")
    nInds = bedfile.shape[0]
    outped = str(args.bfile) + "_pseudoHap.ped"
    for ind_idx in range(nInds):
        print("    working on sample " + str(bedfile.iid.values[ind_idx]) + " (" + str(ind_idx+1) + "/" + str(nInds) + ")")
        make_pseudohap(x=ind_idx, bedfile=bedfile, outped=outped)
    os.system("awk \'{print $1, $2, $3, $4}\' " + str(args.bfile)+ ".bim > " + str(args.bfile) + "_pseudoHap.map")
    args.bfile = str(args.bfile) + "_pseudoHap"
    os.system("plink1.9 --file " + str(args.bfile) + " --make-bed --out " + str(args.bfile))
    os.system("rm " + str(args.bfile) + ".ped")
    os.system("rm " + str(args.bfile)+ ".map")
    return 0

### Define 'make pseudohaploid' function
def make_pseudohap(x, bedfile, outped):
    genos = bedfile.sel(sample=bedfile.sample.values[x]).values.tolist()
    het_idx = np.where(np.array(genos)==1)[0][0:]
    for g in het_idx:
        genos[g] = float(np.array(random.sample([0., 2.], k=1)))
    alleles = [0] * len(genos)
    for n in range(len(genos)):
        if genos[n]==0.:
            alleles[n] = bedfile.a0.values[n]
        elif genos[n]==2.:
            alleles[n] = bedfile.a1.values[n]
    IND_ped = bedfile.fid.values[x] + " " + bedfile.iid.values[x] + " " + bedfile.father.values[x] + " " + bedfile.mother.values[x] + " " + bedfile.gender.values[x] + " " + bedfile.trait.values[x] + str([x for pair in zip(alleles, alleles) for x in pair]).replace("["," ").replace("'","").replace(",","").replace("]","")
    if x==0:
        with open(outped, "w") as pedfile:
            pedfile.write(IND_ped)
    else:
        with open(outped, "a") as pedfile:
            pedfile.write("\n" + IND_ped)

### Run script
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

