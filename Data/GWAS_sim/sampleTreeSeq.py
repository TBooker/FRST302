## A script for grabbing VCF files and phenotypes from Tree-seq files.
# It will a) add neutral mutations
#         b) sample individuals at random from the metapop
#         c) record the phenotypes and location info of individuals
import argparse
import tskit, pyslim, msprime
import numpy as np
import pandas as pd
from localAdaptationLib import *

def round_down(num, divisor):
    return num - (num%divisor)

def main():
## Define command line args
    parser = argparse.ArgumentParser(description="A script to generate VCF files for individuals sampled from simulated populations")

    parser.add_argument("--tree", "-t",
    required = True,
    dest = "tree",
    type = str,
    help = "The file containing the tree sequence of the simulations")

    parser.add_argument("--individuals", "-i",
    required = True,
    dest = "individuals",
    type = int,
    help = "The number of individuals you would like to sample from each population")

    parser.add_argument("--output", "-o",
    required = True,
    dest = "output",
    type = str,
    help = "The output file name - don't give an extension")


    args = parser.parse_args()

    tree_file = args.tree
    file_dir = "./"
    
# Load in the tree sequence
    orig_ts = pyslim.update( tskit.load(args.tree) )
    print([i for i in orig_ts.populations()])
    p1_alive = pyslim.individuals_alive_at(orig_ts, 0, population=1)
    p2_alive = pyslim.individuals_alive_at(orig_ts, 0, population=2)
    print(f"There are {len(p1_alive)} p1 individuals alive in the final generation.")
    print(f"There are {len(p2_alive)} p2 individuals alive in the final generation.")
    ts = orig_ts.simplify()

    effectFile = open(args.output+".muts.csv", "w")
    effectFile.write("id,position,effect\n")
    effectsList = []
    for variant in ts.variants():
        metaDat = ts.site(variant.site.id)
        position = metaDat.position
        effect = metaDat.mutations[0].metadata["mutation_list"][0]["selection_coeff"]
        effectsList.append(effect)
        effectFile.write(",".join(map(str,[variant.site.id, position, effect]))+"\n")
    effectFile.close()
    effectsArray = np.array(effectsList)

    additive = (ts.genotype_matrix().T * effectsArray).sum(axis = 1).reshape(-1,2).sum(axis=-1)
    h_2 = 0.4
    V_A = np.std(additive)**2
    V_E = (V_A - h_2 * V_A) / h_2 #    // from h2 == V_A / (V_A + V_E)
    env = np.random.normal(0.0,  np.sqrt(V_E), len(additive));

    phens = additive+env
    phenOut = open(args.output+".phen.csv", "w")
    phenOut.write("id,phen,pop\n")

    for p in range(1000):
        if p<500:
            phenOut.write(",".join(["tsk_"+str(p),str(phens[p]),"p1"])+"\n")
        if p>=500:
            phenOut.write(",".join(["tsk_"+str(p),str(phens[p]),"p2"])+"\n")
    phenOut.close()



# Mutate the tree sequence - adding neutral polymorphism
    mutated_ts = msprime.sim_mutations(ts,
                rate=3e-8,
                model=msprime.SLiMMutationModel(type=0),
                random_seed=54321,
                keep = True)
    outVCF = open(args.output+".vcf","w")
    mutated_ts.write_vcf(outVCF, 
                            contig_id="chr1")
    outVCF.close()                            
    return

if __name__ == "__main__":
    main()
