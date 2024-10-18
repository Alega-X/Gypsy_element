### tblastn search for GAG-POL elements --- START ---
DIRECTORY="/g/data/lf10/references/Arthropoda_genome_assemblies_2022"
DIRECTORY_scratch="/scratch/lf10/rh1772/Arthropoda_genome_assemblies_2022"
DIRECTORY_Zeyu="/scratch/lf10/zx6715/analysis/tBLASTn"
export PATH="/g/data/lf10/tools/ncbi-blast-2.9.0+-src/rmblast/bin/:${PATH}"
export PATH="/g/data/lf10/tools/EMBOSS-6.6.0/emboss/:${PATH}"
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"

### running tBlastn
#start from Diptera class, analyze 3 classes at a time.
#loop 1: Diptera Lepidoptera Coleoptera
#loop 2: Hymenoptera Hemiptera Polyneoptera
#loop 3: Palaeoptera
#CLASS="Diptera"
#CLASS="Lepidoptera"
#CLASS="Coleoptera"
#CLASS="Hymenoptera"
#CLASS="Hemiptera"
#CLASS="Polyneoptera"
#CLASS="Palaeoptera"

### run a for loop or run individual jobs per ASSEMBLY
for CLASS in Hymenoptera Hemiptera Polyneoptera Palaeoptera; do
cat ${DIRECTORY}/${CLASS}/${CLASS}_reference-genomes.assembly.species.family.size.txt | while read ASSEMBLY SPECIES FAMILY SIZE; do
### common part
qsub -v DIRECTORY=${DIRECTORY},DIRECTORY_scratch=${DIRECTORY_scratch},DIRECTORY_Zeyu=${DIRECTORY_Zeyu},CLASS=${CLASS},ASSEMBLY=${ASSEMBLY},SPECIES=${SPECIES} -P lf10 -l walltime=6:00:00,mem=24GB,ncpus=2,storage=scratch/lf10+gdata/lf10,jobfs=20GB -q expresssr -o /home/150/zx6715/logs/ -e /home/150/zx6715/logs/ /scratch/lf10/zx6715/scripts/GAG-POL_tblastn_ZX1.sh
fi
done
### tblastn search for GAG-POL elements --- END ---


#Job content (GAG-POL_tblastn.ZX1.sh)
export PATH="/g/data/lf10/tools/ncbi-blast-2.9.0+-src/rmblast/bin/:${PATH}"
export PATH="/g/data/lf10/tools/EMBOSS-6.6.0/emboss/:${PATH}"
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"
module load bedtools/2.28.0
mkdir -p ${DIRECTORY_Zeyu}/tblastn_results
mkdir -p ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results
mkdir -p ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/ENV
mkdir -p ${DIRECTORY_Zeyu}/${CLASS}/FNA

### run this common part per ASSEMBLY --- START ---
### blastn database made per genome assembly: ${DIRECTORY}/${CLASS}/blastdb/${ASSEMBLY}_NCBI_blast
for ORF in GAG POL; do
fasta_file="/scratch/lf10/zx6715/Data/baits_sum_${ORF}.fasta"
tblastn -db ${DIRECTORY}/${CLASS}/blastdb/${ASSEMBLY}_NCBI_blast -outfmt 10 -num_threads 6 \
-query ${fasta_file} \
-out ${DIRECTORY_Zeyu}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_tblastn_results.out.csv
done

### summarising the 1st and 2nd tBlastn results --- START ---
for ORF in GAG POL; do
### extracting high-score GAG POL insertions and merge entries
### longer than 300nt for GAG and longer than 1500nt for POL
if [[ ${ORF} == "POL" ]]; then
### extract ORFs tblastn hits
cat ${DIRECTORY_Zeyu}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_tblastn_results.out.csv | tr ',' ' ' |\
awk '{if ($9<$10 && $12>50) print $2,$9-1,$10,$1":"$7"-"$8,$5":"$6":"$11":"$12,"+";
else if ($9>$10 && $12>50) print $2,$10-1,$9,$1":"$7"-"$8,$5":"$6":"$11":"$12,"-"
}' | tr ' ' '\t' | sort -k1,1 -k2,2n > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.bed

bedtools merge -d 300 -s -c 5,6 -o count,distinct -i ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.bed |\
awk '{if ($3-$2>1500 && $4>1) print $1,$2,$3,$4,".",$5}' | tr ' ' '\t' > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.high-score.bed
elif [[ ${ORF} == "GAG" ]]; then
### extract ORFs tblastn hits
cat ${DIRECTORY_Zeyu}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_tblastn_results.out.csv | tr ',' ' ' |\
awk '{if ($9<$10 && $12>30) print $2,$9-1,$10,$1":"$7"-"$8,$5":"$6":"$11":"$12,"+";
else if ($9>$10 && $12>30) print $2,$10-1,$9,$1":"$7"-"$8,$5":"$6":"$11":"$12,"-"
}' | tr ' ' '\t' | sort -k1,1 -k2,2n > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.bed

bedtools merge -d 300 -s -c 5,6 -o count,distinct -i ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.bed |\
awk '{if ($3-$2>300 && $4>1) print $1,$2,$3,$4,".",$5}' | tr ' ' '\t' > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.high-score.bed
fi
done

### obtain coordinate of high score ENV insertions, extended by 2000nt at either end
awk '{print $1,$2-2000,$3+2000,$4,$5,$6}' ${DIRECTORY_scratch}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_1st-and-2nd-and-3rd-and-4th-and-5th_ENV_preset_tblastn_results.out.high-score.bed |\
awk '{if($2<0) $2=0; print}' | tr ' ' '\t' > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/ENV/${ASSEMBLY}_${SPECIES}_1st-and-2nd-and-3rd-and-4th-and-5th_ENV_preset_tblastn_results.out.high-score.2000nt_extended.bed

### collect GAG-POL insertions that do not have ENV next to them.
awk '{if($6=="+") print $1,$2,$3+500,$4,$5,$6; else print $1,$2-500,$3,$4,$5,$6}' ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG_preset_tblastn_results.out.high-score.bed |\
awk '{if($2<0) $2=0; print}' | tr ' ' '\t' |\
bedtools intersect -s -wao -a - -b ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_POL_preset_tblastn_results.out.high-score.bed |\
awk '{if($6=="+" && $NF!=0) print $1,$2,$9,$4":"$10,".",$6; else if($6=="-" && $NF!=0) print $1,$8,$3,$4":"$10,".",$6}' |\
awk '{if($3-$2>3000 && $3-$2<10000) print $1,$2,$3,$1":"$2":"$3":"$4,".",$6}' | tr ' ' '\t' |\
bedtools intersect -v -a - -b ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/ENV/${ASSEMBLY}_${SPECIES}_1st-and-2nd-and-3rd-and-4th-and-5th_ENV_preset_tblastn_results.out.high-score.2000nt_extended.bed > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.bed

#check whether corresponding .fna file exists.
if [[ ! -f "${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna" ]]; then
zcat ${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna.gz > ${DIRECTORY_Zeyu}/${CLASS}/FNA/${ASSEMBLY}_genomic.fna
Fna_File="${DIRECTORY_Zeyu}/${CLASS}/FNA/${ASSEMBLY}_genomic.fna"
else
Fna_File="${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna"
fi

for ORF in GAG POL; do
bedtools intersect -wo -s -a ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.bed \
-b ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.high-score.bed |\
awk -v ORF=${ORF} '{print $7,$8,$9,$4":"ORF,".",$12}' | tr ' ' '\t' |\
bedtools getfasta -name -s -fi ${Fna_File} -bed - | fasta_formatter - -t | tr ':' '@' |\
awk '{print ">"$1"\n"toupper($2)}' > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.${ORF}.fasta

transeq ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.${ORF}.fasta ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.${ORF}.pep2
cat ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.${ORF}.pep2 | tr '@' ':' | fasta_formatter - -t |\
awk '{if (gsub("*","*",$2)<5) print $1,$2,length($2)}' | sort -k3,3nr | awk '!seen[$1]++' | awk '{print ">"$1"\n"$2}' > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.${ORF}.pep
done

### only take GAG and POL that both survived in the same locus
cat ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.*.pep |\
fasta_formatter - -t | awk '{split($1,a,":"); GYPSY[a[1]":"a[2]":"a[3]]++} END {for(var in GYPSY) {if(GYPSY[var]==2) print var}
}' > ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.GAG-and-POL.list

for ORF in GAG POL; do
cat ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.GAG-and-POL.list | while read NAME; do
fasta_formatter -i ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.${ORF}.pep -t |\
awk -v NAME=${NAME} '{if($1~NAME) print ">"$1"\n"$2}' >> ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.${ORF}.GAG-and-POL.pep
done
done
### run this common part per ASSEMBLY --- END ---
