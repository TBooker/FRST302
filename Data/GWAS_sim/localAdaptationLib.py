## A library of Python functions that I wrote for the LocalAdaptationArchitechture repo
import numpy as np

class Mutation:
  def __init__(self, raw, effect_dictionary):
    self.ID = int( raw[0] )
    self.permID = int( raw[1] )
    self.mutType = raw[2]
    self.position = int(raw[3])
    self.effectSize = float(raw[4])
    self.dominance = float(raw[5])
    self.sourcePop = raw[6]
    self.generation = int(raw[7])
    self.numCopies = int(raw[8])
#    print(raw)
#    print(raw[4][0:6])
#    print(effect_dictionary)
    self.effect_1 = effect_dictionary[self.permID]["e1"]
    self.effect_2 = effect_dictionary[self.permID]["e2"]

class Individual:
    def __init__(self, raw, genome_dict, optima, mutDCT, maskMutationIndex):

        self.rawID = raw[0]
        self.pop = raw[0].split(":")[0]
        self.indID = raw[0].split(":")[1]
        self.sex =  raw[1]
        self.genome_1_raw = raw[2]
        self.genome_2_raw = raw[3]

        self.genome_1 = genome_dict[self.genome_1_raw]
        self.genome_2 = genome_dict[self.genome_2_raw]

        self.genome = np.array([self.genome_1.mutArray, self.genome_2.mutArray])

## Actually read in the list of environments though
        self.environmentalOptimum = optima[ self.pop ]
        if maskMutationIndex==-1:
            self.phenotype = self.calcPhenotype(mask = maskMutationIndex,)
        else:
# This is the panmictic allele frequency
            maskMutFreq = mutDCT[maskMutationIndex].numCopies/(196*100*2)
            maskMutEffect = mutDCT[maskMutationIndex].effect_1
#            print("\n\n$$$", maskMutationIndex,maskMutFreq, maskMutEffect)
# Calculate the phenotype
            self.phenotype = self.calcPhenotype(mask = maskMutationIndex, alleleFreq=maskMutFreq, mut_effect = maskMutEffect)

# This method converts the genotype at the masking locus to homozygous wildtype
    def calcPhenotype_zero(self, mask=-1):
        assert type(mask)==int
        assert mask >=-1

        if mask == -1:
            return self.genome.sum(axis = 0).sum()
        else:
            zero_mask = np.zeros(self.genome.shape[1])
            zero_mask[mask] = 1
            masker = zero_mask==0

            return self.genome.transpose()[masker].transpose().sum(axis = 0).sum()
# This method converts the genotype at masking locus to one based on random raws of alleles
    def calcPhenotype(self, mask=-1,alleleFreq=0, mut_effect = 0):
        assert type(mask)==int
        assert mask >=-1
# If not masking the mutation, just return the actual phenotype
        if mask == -1:
            return self.genome.sum(axis = 0).sum()
        else:
# Replace the mutation entry in the genome array with a effects that are randomly added
# Based on the panmictic allele freq.
            new_genotype = mut_effect* (np.random.uniform(size =2 )<=alleleFreq)
            new_genome = self.genome.copy().transpose()
            new_genome[mask] = new_genotype
            return new_genome.transpose().sum(axis = 0).sum()
    def calcFitness(self, optimum, mask = -1):
        if mask == -1:
            phenotype = self.calcPhenotype()
        else:
            phenotype = self.calcPhenotype(mask=mask)
#        print(self.rawID)
#        print(phenotype)

        part_1 = ((optimum - phenotype)**2)/(2*15.)
#        print(optimum - phenotype)
#        print(part_1)
        relFitness = np.exp(-1*part_1);
#        print(relFitness)
        return relFitness

class Genome:
  def __init__(self, raw, mutationDict):
    self.rawID = raw[0]
    self.pop = raw[0].split(":")[0]
    self.ID = int( raw[0].split(":")[1] )
    self.chromType =  raw[1]
    self.mutList_raw = [int(k) for k in  raw[2:] ]
    full_mutation_list = np.zeros(len(  mutationDict.keys() ))
    for i in self.mutList_raw:
        full_mutation_list[i] = mutationDict[i].effect_1
    self.mutArray = full_mutation_list
#    print("!",self.mutList_raw)

def grab_lines(input_file_, start, stop ):
    shall_I_yield = False
    line_container = []
    for l in open( input_file_, "r"):

        line = l.strip()
        if line.startswith(stop):
            shall_I_yield = False
        if shall_I_yield == False:
            pass
        elif shall_I_yield == True:
            line_container.append( line.split(" ") )

        if line.startswith(start):
            shall_I_yield = True
    if len(line_container) == 0:
        return None
    else:
        return line_container

def grab_mutations(input_file, effectDict):
    mutation_list = [Mutation(m, effectDict) for m in grab_lines(input_file, "Mutations:", "Individuals:" ) ]
    mutation_dict = {m.ID: m for m in mutation_list}

    return( mutation_dict )

def grab_individuals(input_file, genome_dictionary, optimaDict, mutation_dictionary,  maskingMutation):
    individuals_list = [Individual(ind, genome_dictionary, optimaDict, mutation_dictionary, maskingMutation) for ind in grab_lines(input_file, "Individuals:", "Genomes:" )]
    return( individuals_list)

def grab_genomes(input_file, mutations):
    genome_list = [Genome(g, mutations) for g in grab_lines(input_file, "Genomes:", "COOLSTUFF" ) ]
    genome_dict = {g.rawID: g for g in genome_list}
    return( genome_dict )

def get_effect_dict(effect_file_raw):
    effectDict = {}
    for e_raw in open(effect_file_raw):
        if e_raw.startswith("ID,position,selCoeff"):continue
        e = e_raw.strip().split(",")

        effectDict[int(float(e[0]))] = {"effect":float(e[2]),
        "position":int(float(e[1])),
        "e1":float(e[3]),
        "e2":float(e[4])}
    return(effectDict)


def down_sample(x, f=7):
    # pad to a multiple of f, so we can reshape
    # use nan for padding, so we needn't worry about denominator in
    # last chunk
    xp = np.r_[x, np.nan + np.zeros((-len(x) % f,))]
    # reshape, so each chunk gets its own row, and then take mean
    return np.nanmean(xp.reshape(-1, f), axis=-1)
