#!/bin/bash
for i in {1..688}; do
ORF=$(awk -v i="$i" 'NR==i {gsub(/>/, ""); print}' /scratch/lf10/zx6715/Data/POL_real.txt)
hhr_file="/scratch/lf10/zx6715/analysis/hhr/${ORF}.hhr"
a3m_file="/scratch/lf10/zx6715/analysis/a3m/${ORF}.a3m"
input_fasta="/scratch/lf10/zx6715/HHdata/${ORF}.fasta"

if [[ ! -f ${hhr_file} ]]; then
qsub -v hhr_file=${hhr_file},a3m_file=${a3m_file},input_fasta=${input_fasta} -P lf10 -l walltime=6:00:00,mem=24GB,ncpus=2,storage=scratch/lf10+gdata/lf10,jobfs=20GB -q expresssr -o /home/150/zx6715/logs/ -e /home/150/zx6715/logs/ /scratch/lf10/zx6715/scripts/hhZX.sh
fi
done

#The job content (hhZX.sh)
### first we need to make a3m file
/g/data/lf10/tools/hh-suite/bin/hhblits -i ${input_fasta} -d /g/data/lf10/tools/hh-suite/databases/UniRef30_2022_02 -oa3m ${a3m_file} -n 1

### same parameter as HHpred web search
/g/data/lf10/tools/hh-suite/bin/hhsearch -i ${a3m_file} -o ${hhr_file} -d /g/data/lf10/tools/hh-suite/pdb70 \
-p 20 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 \
-contxt /g/data/lf10/tools/hh-suite/data/context_data.crf
