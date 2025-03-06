cd GWAS
## Assuming that you've run SLiM and have the tree-seq output
python ../sampleTreeSeq.py --tree ../twoPopModelForTeachingGWAS.X.ts --individuals 10 --output gwasDemo 

# Edit VCF file from msprime/tskit
python ../editVCF.py gwasDemo.vcf gwasDemo.edit.vcf
mv gwasDemo.edit.vcf gwasDemo.vcf

 
# Generate plink files from VCF
~/software/plink_mac/plink --make-bed --vcf gwasDemo.vcf --out gwasDemo --set-missing-var-ids @:# --double-id --allow-extra-chr

# Add phenotypes to the .fam file 
python3 ../editPlinkFamFile.py gwasDemo.fam twoPopModelForTeachingGWAS.phen.txt gwasDemo.phen.tweak.fam
mv gwasDemo.phen.tweak.fam gwasDemo.fam  

# Run GWAS on all SNPs at once...
# Generate kinship matrix
~/software/gemma \
 -bfile gwasDemo \
 -gk \
 -o gwasDemo


 ~/software/gemma \
 -bfile gwasDemo \
 -k output/gwasDemo.cXX.txt  \
 -lmm 4 \
 -o gwasDemo 
