import sys
import pandas as pd

# Read the plink file
fam = pd.read_csv(sys.argv[1], header = None, sep='\s+')

# Read the phenotype file
phen = pd.read_csv(sys.argv[2], header = None, names=["phen"])
phen["id"] = ["tsk_"+str(i) for i in range(len(phen.phen))]

phen_dict =  pd.Series(phen.phen.values,index=phen.id).to_dict()

#Map the phenotypes by individual ID
fam[5] = fam[1].map(phen_dict)

# Output the file
fam.to_csv(sys.argv[3],
			header = None,
			sep = " ",
			index = False)
