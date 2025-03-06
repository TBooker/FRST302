# A shell script to prepare the tskit VCFs for a GWAS...

vcf=$1
phen=$2
run=$3
outPref=$4
vcf_func=$5

plink=~/software/plink_mac/plink
$plink --make-bed \
	--vcf $vcf_neu \
	--out $outPref/nucl.$run \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr

	$plink \
		--bfile $outPref/nucl.$run \
		--indep-pairwise 100 10 0.2 \
		--out $outPref/nucl.$run \
		--make-founders \
		--allow-extra-chr

		$plink \
		--bfile $outPref/nucl.$run \
		--extract $outPref/nucl.${run}.prune.in \
		--make-bed \
		--out $outPref/nucl.${run}.LDpruned \
		--allow-extra-chr

script_bin=~/UBC/LocalAdaptationArchitechture/bin/
# Use a little python script to edit the fam file to include the phenotypes...
python3 $script_bin/editPlinkFamFile.py $outPref/nucl.${run}.LDpruned.fam $phen $outPref/nucl.${run}.LDpruned.tweak.fam 15
mv $outPref/nucl.${run}.LDpruned.tweak.fam $outPref/nucl.${run}.LDpruned.fam

# Run the GWAS
	#first use the option -gk to generate the relatedness matrix
gemma \
 -bfile $outPref/nucl.$3.LDpruned \
 -gk \
 -o nucl.$3.LDpruned \
 -outdir $outPref/output

# Generate Plink files for the functional sites:

$plink --make-bed \
	--vcf $vcf_func \
	--out $outPref/nuclFun.$run \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr

python3 $script_bin/editPlinkFamFile.py $outPref/nuclFun.${run}.fam $phen $outPref/nuclFun.${run}.tweak.fam 15
mv $outPref/nuclFun.${run}.tweak.fam $outPref/nuclFun.${run}.fam


#run the GWA controlling for relatedness, with the output matrix (*.cXX.txt) only on those sites that are functional
gemma \
 -bfile $outPref/nuclFun.$run \
 -k $outPref/output/nucl.${run}.LDpruned.cXX.txt \
 -lmm 4 \
 -o nucl.${run}.LDpruned \
 -outdir $outPref/output
