#!/bin/bash
# Use: bash assemblyCASPER.sh NOMBRE_TRABAJO

SEQS=$(pwd)
SALIDAS=$(pwd)
BIN=/usr/local/bin
BIN2=/usr/local/bin
BIN3=/usr/bin
COUNT=0
for FAA in `cat lista`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr

echo "$BIN/casper <(zcat $SEQS/$FAA"_R1.fastq.gz") <(zcat $SEQS/$FAA"_R2.fastq.gz") -o $FAA.assembly.fastq -o $FAA"_assembly"" >>$*.$COUNT.scr

echo "$BIN3/fastq_to_fasta -n -i $FAA.assembly.fastq -o$FAA"_assembly.fas"" >>$*.$COUNT.scr
done
