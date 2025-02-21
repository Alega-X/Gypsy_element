#!/bin/bash
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"

DIRECTORY="/g/data/lf10/references/Arthropoda_genome_assemblies_2022"
DIRECTORY_Zeyu="/scratch/lf10/zx6715/analysis"
MMSEQS="/scratch/lf10/zx6715/analysis/mmseqs"

#We start from Palaeoptera, then: Polyneoptera, Hemiptera, Hymenoptera, Coleoptera, Lepidoptera
#There is too many elements in Lepidoptera that the qsub command cannot analyze all elements at a time. We need to submit again for elements left.
#CLASS="Diptera"
CLASS="Lepidoptera"
#CLASS="Coleoptera"
#CLASS="Hymenoptera"
#CLASS="Hemiptera"
#CLASS="Polyneoptera"
#CLASS="Palaeoptera"

MMRESULT="${MMSEQS}/${CLASS}"
mkdir -p ${DIRECTORY_Zeyu}/rep_hhpred/hhr/${CLASS}

#Read the simplified assemblies list to go over all of species for each class
cat ${DIRECTORY_Zeyu}/Class_list/${CLASS}_list.txt | while read ASSEMBLY SPECIES FAMILY SIZE; do

#For species containing representative elements, read the summary file, extract sequence for each single element and submit a HHpred job
if [[ -f "${DIRECTORY_Zeyu}/cluster_rep/${CLASS}/${SPECIES}.fasta" ]]; then
NUM=$(($(wc -l < ${DIRECTORY_Zeyu}/cluster_rep/${CLASS}/${SPECIES}.fasta) / 2))
for N in $(seq 1 $NUM); do
NAME=`(fasta_formatter -i ${DIRECTORY_Zeyu}/cluster_rep/${CLASS}/${SPECIES}.fasta -t |\
grep "${SPECIES}_${N}:POL" | awk '{print $1}')`

#Check whether the element have been extracted, if not, generate a fasta file for this element
if [[ ! -f ${DIRECTORY_Zeyu}/rep_hhpred/hhr/${CLASS}/${NAME}.fasta ]]; then
SEQUENCE=`(fasta_formatter -i ${DIRECTORY_Zeyu}/cluster_rep/${CLASS}/${SPECIES}.fasta -t |\
grep "${SPECIES}_${N}:POL" | awk '{print $2}')`
echo ">"${NAME} >> ${DIRECTORY_Zeyu}/rep_hhpred/hhr/${CLASS}/${NAME}.fasta
echo ${SEQUENCE} >> ${DIRECTORY_Zeyu}/rep_hhpred/hhr/${CLASS}/${NAME}.fasta
fi

#Set up for HHpred
hhr_file="${DIRECTORY_Zeyu}/rep_hhpred/hhr/${CLASS}/${NAME}.hhr"
a3m_file="${DIRECTORY_Zeyu}/rep_hhpred/hhr/${CLASS}/${NAME}.a3m"
input_fasta="${DIRECTORY_Zeyu}/rep_hhpred/hhr/${CLASS}/${NAME}.fasta"
if [[ ! -f ${hhr_file} ]]; then
qsub -v hhr_file=${hhr_file},a3m_file=${a3m_file},input_fasta=${input_fasta} -P lf10 -l walltime=6:00:00,mem=24GB,ncpus=2,storage=scratch/lf10+gdata/lf10,jobfs=20GB -q expresssr -o /home/150/zx6715/logs/ -e /home/150/zx6715/logs/ /scratch/lf10/zx6715/scripts/hhZX.sh
fi

done
fi
done

#Job content (hhZX.sh)
### first we need to make a3m file
/g/data/lf10/tools/hh-suite/bin/hhblits -i ${input_fasta} -d /g/data/lf10/tools/hh-suite/databases/UniRef30_2022_02 -oa3m ${a3m_file} -n 1

### same parameter as HHpred web search
/g/data/lf10/tools/hh-suite/bin/hhsearch -i ${a3m_file} -o ${hhr_file} -d /g/data/lf10/tools/hh-suite/pdb70 \
-p 20 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 \
-contxt /g/data/lf10/tools/hh-suite/data/context_data.crf
