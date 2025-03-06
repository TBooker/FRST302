# Edit the VCF from the simulation
import sys

vcf = open(sys.argv[1])

snp_counter = 0

new_vcf=open(sys.argv[2],"w")

for l in vcf:
    if l.startswith("#"):
        new_vcf.write(l)
        continue
  #  else:break

    snp_counter+=1
    v = l.split("\t")
    v[2]="snp"+str(snp_counter)
    v[3]="A"
    v[4]="T"
    new_vcf.write("\t".join(v)+"\n")
new_vcf.close()