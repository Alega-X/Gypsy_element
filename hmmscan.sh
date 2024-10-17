### this code is to identify GAG-POL elements in dipteran genomes --- 2024-Aug-05 ---

DIRECTORY="/g/data/lf10/references/Arthropoda_genome_assemblies_2022"
DIRECTORY_scratch="/scratch/lf10/rh1772/Arthropoda_genome_assemblies_2022"

mkdir -p ${DIRECTORY}/GAG-POL_Diptera
export PATH="/g/data/lf10/tools/fastx_toolkit-0.0.14/bin/:${PATH}"

### extract CDS from embl file
awk '{if($1=="FT" && NF==2) print $NF}' ${DIRECTORY}/GAG-POL_Diptera/Drosophila_gypsy_RepBase_2024-Aug.embl |\
grep -v ")" | grep -v "," | tr -d '=/"' | sed 's/translation//g' | sed 's/notegag.//g' | sed 's/noteGAG.//g' | sed 's/noteGag.//g' | sed 's/noteenv.//g' | sed 's/noteENV.//g' | sed 's/notepol.//g' | sed 's/notePOL.//g' | sed 's/product/>/g' |\
fasta_formatter - -t | awk '{print ">"$1"\n"$2}' > ${DIRECTORY}/GAG-POL_Diptera/Drosophila_gypsy_RepBase_2024-Aug.CDS.fasta

### hmmscan
export PATH="/g/data/lf10/tools/hmmer-3.3.1/bin/:${PATH}"
mkdir -p ${DIRECTORY}/GAG-POL_Diptera/hmmscan
mkdir -p ${DIRECTORY}/GAG-POL_Diptera/HHpred
fasta_formatter -i ${DIRECTORY}/GAG-POL_Diptera/Drosophila_gypsy_RepBase_2024-Aug.CDS.fasta -t | while read gene AA; do
fasta_formatter -i ${DIRECTORY}/GAG-POL_Diptera/Drosophila_gypsy_RepBase_2024-Aug.CDS.fasta -t |\
awk -v gene=${gene} '{if($1==gene) print ">"gene"\n"$2}' > ${DIRECTORY}/GAG-POL_Diptera/HHpred/${gene}.fasta
hmmscan -o ${DIRECTORY}/GAG-POL_Diptera/hmmscan/${gene}.hhr /g/data/lf10/tools/hh-suite/GyDB_collection/database/GyDB ${DIRECTORY}/GAG-POL_Diptera/HHpred/${gene}.fasta
awk -v gene=${gene} '{if(NR>15 && NR<26) print gene,$NF,NR-15}' ${DIRECTORY}/GAG-POL_Diptera/hmmscan/${gene}.hhr >> ${DIRECTORY}/GAG-POL_Diptera/hmmscan/hmmscan.Drosophila_gypsy_RepBase_2024-Aug.stats
done
