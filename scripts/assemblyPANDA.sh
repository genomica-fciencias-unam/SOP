#!/bin/bash
# Use: bash assembly.sh NOMBRE_TRABAJO

SEQS=$(pwd)
SALIDAS=$(pwd)
BIN=/usr/bin
BIN2=/usr/bin
COUNT=0

for FAA in `cat lista`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "zcat $SEQS/$FAA"_R1.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr_R1.fastq" &" >>$*.$COUNT.scr 
echo "zcat $SEQS/$FAA"_R2.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr_R2.fastq"" >>$*.$COUNT.scr
echo "$BIN/pandaseq -B -f $SEQS/$FAA"tr_R1.fastq" -r $SEQS/$FAA"tr_R2.fastq" -t 0.95 -l 250 -L 470 -o 15 -w $SALIDAS/$FAA"_$*.fasta" -G $SALIDAS/$FAA"-$*.log.bz2"" >>$*.$COUNT.scr

done
