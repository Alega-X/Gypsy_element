### tblastn search for GAG-POL elements --- START ---
DIRECTORY="/g/data/lf10/references/Arthropoda_genome_assemblies_2022"
DIRECTORY_scratch="/scratch/lf10/rh1772/Arthropoda_genome_assemblies_2022"
DIRECTORY_Zeyu="/scratch/lf10/zx6715/analysis/tBLASTn"
export PATH="/g/data/lf10/tools/ncbi-blast-2.9.0+-src/rmblast/bin/:${PATH}"
export PATH="/g/data/lf10/tools/EMBOSS-6.6.0/emboss/:${PATH}"
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"

### running tBlastn
#start from Diptera class, we do 3 classes a time. some of fna file is missed in first round. here we process the following works
#loop 1: Diptera Lepidoptera Coleoptera round2 for fna
#loop 2: Hymenoptera Hemiptera Polyneoptera round2 for fna + Palaeoptera
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
if [[ ! -f "${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna" ]]; then
qsub -v DIRECTORY=${DIRECTORY},DIRECTORY_scratch=${DIRECTORY_scratch},DIRECTORY_Zeyu=${DIRECTORY_Zeyu},CLASS=${CLASS},ASSEMBLY=${ASSEMBLY},SPECIES=${SPECIES} -P lf10 -l walltime=6:00:00,mem=24GB,ncpus=2,storage=scratch/lf10+gdata/lf10,jobfs=20GB -q expresssr -o /home/150/zx6715/logs/ -e /home/150/zx6715/logs/ /scratch/lf10/zx6715/scripts/GAG-POL_tblastn_fna.sh
fi
done
done

### tblastn search for GAG-POL elements --- END ---

#Job content (GAG-POL_tblastn.fna.sh)
export PATH="/g/data/lf10/tools/ncbi-blast-2.9.0+-src/rmblast/bin/:${PATH}"
export PATH="/g/data/lf10/tools/EMBOSS-6.6.0/emboss/:${PATH}"
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"
module load bedtools/2.28.0
mkdir -p ${DIRECTORY_Zeyu}/tblastn_results
mkdir -p ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results
mkdir -p ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/ENV
mkdir -p ${DIRECTORY_Zeyu}/${CLASS}/FNA

if [[ ! -f "${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna" ]]; then
zcat ${DIRECTORY}/${CLASS}/${ASSEMBLY}_genomic.fna.gz > ${DIRECTORY_Zeyu}/${CLASS}/FNA/${ASSEMBLY}_genomic.fna
fi

for ORF in GAG POL; do
bedtools intersect -wo -s -a ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.bed \
-b ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_${ORF}_preset_tblastn_results.out.high-score.bed |\
awk -v ORF=${ORF} '{print $7,$8,$9,$4":"ORF,".",$12}' | tr ' ' '\t' |\
bedtools getfasta -name -s -fi ${DIRECTORY_Zeyu}/${CLASS}/FNA/${ASSEMBLY}_genomic.fna -bed - | fasta_formatter - -t | tr ':' '@' |\
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
