DIRECTORY="/g/data/lf10/references/Arthropoda_genome_assemblies_2022"
DIRECTORY_Zeyu="/scratch/lf10/zx6715/analysis/tBLASTn"
MMSEQS="/scratch/lf10/zx6715/analysis/mmseqs"

#Analyze 3 or less classes at a time, we start from Palaeoptera
#Second round, Hemiptera Polyneoptera
#Third round, run Coleoptera Hymenoptera
#Final round, run Diptera Lepidoptera
#CLASS="Diptera"
#CLASS="Lepidoptera"
#CLASS="Coleoptera"
#CLASS="Hymenoptera"
#CLASS="Hemiptera"
#CLASS="Polyneoptera"
#CLASS="Palaeoptera"

for CLASS in Diptera Lepidoptera; do
cat ${DIRECTORY}/${CLASS}/${CLASS}_reference-genomes.assembly.species.family.size.txt | while read ASSEMBLY SPECIES FAMILY SIZE; do
if [[ -f "${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.POL.GAG-and-POL.pep" ]]; then
qsub -v DIRECTORY=${DIRECTORY},DIRECTORY_Zeyu=${DIRECTORY_Zeyu},CLASS=${CLASS},ASSEMBLY=${ASSEMBLY},SPECIES=${SPECIES},MMSEQS=${MMSEQS} -P lf10 -l walltime=6:00:00,mem=24GB,ncpus=2,storage=scratch/lf10+gdata/lf10,jobfs=20GB -q expresssr -o /home/150/zx6715/logs/ -e /home/150/zx6715/logs/ /scratch/lf10/zx6715/scripts/cluster_job.sh
fi
done
done

#Job content -start- (cluster_job.sh)
export PATH="/home/150/zx6715/miniconda3/bin:$PATH"
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"

mkdir -p ${MMSEQS}/${CLASS}
MMRESULT="${MMSEQS}/${CLASS}"

### concatenate GAG and POL pep sequeences into one pep file
fasta_formatter -i ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.POL.GAG-and-POL.pep -t |\
awk '{split($1,a,":"); if(a[6] ~ "-") print a[1]":"a[2]":"a[3]":"a[4]":"a[5]" -";
else print a[1]":"a[2]":"a[3]":"a[4]":"a[5]" +"}' | while read COORDINATE STRAND; do
GAG=`(fasta_formatter -i ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.GAG.GAG-and-POL.pep -t |\
grep ${COORDINATE} | awk '{print $2}')`
POL=`(fasta_formatter -i ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.POL.GAG-and-POL.pep -t |\
grep ${COORDINATE} | awk '{print $2}')`
GAG_POL="${GAG} ${POL}"
echo ">"${COORDINATE}":"${STRAND}":GAG_POL" >> ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.GAG_POL.GAG-and-POL.pep
echo ${GAG_POL} >> ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.GAG_POL.GAG-and-POL.pep
done
mmseqs easy-cluster ${DIRECTORY_Zeyu}/${CLASS}/tblastn_results/${ASSEMBLY}_${SPECIES}_GAG-POL_preset.no-nearby-ENV.GAG_POL.GAG-and-POL.pep ${MMRESULT}/${ASSEMBLY}_${SPECIES} /home/150/zx6715/tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1
#Job content -end-
