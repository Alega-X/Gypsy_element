input_fasta="/scratch/lf10/zx6715/analysis/diptera/diptera.RT_extended.commandline_hhpred.fasta"
output_mafft_fasta="/scratch/lf10/zx6715/analysis/diptera/mafft_output.fasta"

### make a multiple sequence alignment using mafft
export PATH="/g/data/lf10/tools/mafft-7.505/bin/:${PATH}"

mafft --auto --thread 8 ${input_fasta} > ${output_mafft_fasta}

### you can check the alignment on mview
https://www.ebi.ac.uk/jdispatcher/msa/mview?stype=protein

### once you check that, the next step is to trim bases at either end of the alignment to use only the part that aligns well
START_POSITION="starting position of the conserved part"
LENGTH="length of the conserved part"
fasta_formatter -i ${output_mafft_fasta} -t | awk '{print ">"$1"\n"substr($2,START,LENGTH)}' > ${output_mafft_trimmed_fasta}

### perform phylogenetic analysis using iqtree
export PATH="/g/data/lf10/tools/iqtree-1.6.12-Linux/bin/:${PATH}"

iqtree -redo -s ${output_mafft_trimmed_fasta} -m rtREV+R4 -bb 1000 -nt 8
