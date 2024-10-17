### fetch the sequence of reverse transcriptase (RT) from the POL ORF.
### this paper describes the domain structure of RT
###https://pubmed.ncbi.nlm.nih.gov/15130474/
### we try to take the part spanning regions that correspond to their fingers, palm, thumb and connection subdomains
###fetch_RT_from_hhr_ZX

module load bedtools/2.28.0
DIRECTORY="/scratch/lf10/zx6715/analysis"
CLASS="diptera"
DOMAIN="RT_extended"
if [[ -f ${DIRECTORY}/${CLASS}/${CLASS}.${DOMAIN}.commandline_hhpred.fasta ]]; then
rm ${DIRECTORY}/${CLASS}/${CLASS}.${DOMAIN}.commandline_hhpred.fasta
fi
awk '
/^>/ {header=$0; next}
{
    sequence = $0
    print header, sequence, length(sequence)
}' ${DIRECTORY}/${CLASS}/${CLASS}.POL.fasta | while read gene AA SIZE; do
gene=$(echo "$gene" | sed 's/>//')
FILE="${DIRECTORY}/hhr/${gene}.hhr"
if [[ -f $FILE ]]; then
### RT, including hits that cover RT regions of 4MH8_A and 4G10_B allowing them to miss 30aa at either ends
awk -v gene=${gene} '{split($12,a,"-"); split($11,b,"-"); if ($2=="4MH8_A" && NF==13 && a[1]<59 && a[2]>395) print gene,b[1]-a[1]+17,b[2]-a[2]+451;
else if ($2=="4G1Q_B" && NF==13 && a[1]<38 && a[2]>374) print gene,b[1]-a[1]-4,b[2]-a[2]+430}' $FILE |\
awk -v SIZE=${SIZE} '{if($3>SIZE) $3=SIZE; if($2<0) $2=0; print}'
fi
done | tr ' ' '\t' | sort -k1,1 -k2,2n |\
bedtools merge -i - |\
awk '{print $1,$2,$3,$1"_RT . +"}' | tr ' ' '\t' |\
bedtools getfasta -name -fi ${DIRECTORY}/${CLASS}/${CLASS}.POL.fasta -bed - >> ${DIRECTORY}/${CLASS}/${CLASS}.${DOMAIN}.commandline_hhpred.fasta
