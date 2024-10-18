#This is to isolate representative elements from mmseq results
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"

DIRECTORY="/g/data/lf10/references/Arthropoda_genome_assemblies_2022"
DIRECTORY_Zeyu="/scratch/lf10/zx6715/analysis"
MMSEQS="/scratch/lf10/zx6715/analysis/mmseqs"

#Analyze one class at a time, start from Palaeoptera
#Then: Polyneoptera, Hemiptera, Hymenoptera, Coleoptera, Lepidoptera, Diptera
CLASS="Diptera"
#CLASS="Lepidoptera"
#CLASS="Coleoptera"
#CLASS="Hymenoptera"
#CLASS="Hemiptera"
#CLASS="Polyneoptera"
#CLASS="Palaeoptera"

MMRESULT="${MMSEQS}/${CLASS}"
mkdir -p ${DIRECTORY_Zeyu}/cluster_rep/${CLASS}

#Identify representatives that hit >=5 elements
cat ${DIRECTORY}/${CLASS}/${CLASS}_reference-genomes.assembly.species.family.size.txt | while read ASSEMBLY SPECIES FAMILY SIZE; do
if [[ -f "${MMRESULT}/${ASSEMBLY}_${SPECIES}_cluster.tsv" ]]; then
NUM=1
awk '{print $1}' ${MMRESULT}/${ASSEMBLY}_${SPECIES}_cluster.tsv | sort | uniq -c | awk '$1 >= 5 {print $2}'| while read REP; do
SEQ=`(fasta_formatter -i ${MMRESULT}/${ASSEMBLY}_${SPECIES}_rep_seq.fasta -t |\
grep ${REP} | awk '{print $3}')`
echo ">"${SPECIES}"_"${NUM}":POL" >> ${DIRECTORY_Zeyu}/cluster_rep/${CLASS}/${SPECIES}.fasta
echo $SEQ >> ${DIRECTORY_Zeyu}/cluster_rep/${CLASS}/${SPECIES}.fasta
NUM=$((NUM+1))
done
fi
done
