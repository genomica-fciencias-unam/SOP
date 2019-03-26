Para los metagenomas:

Los pasos generales son: 

<<<<<<< HEAD
0. Mapeo de secuencias contra el hospedero (Bowtie2)
```
comando del mal
=======
0. Mapeo de secuencias contra el hospedero utilizando Bowtie2:

a) Se selecciona un genoma de referencia para realizar el filtrado, es necesario formatearlo:
```
/srv/home/mromero/mbin/bowtie2/bowtie2-build /srv/home/chernandez/processing/paired/GCA_002806865.2_ASM280686v2_genomic.fna calabacitagenome
```
b) Realizar el alineamiento contra el genoma de referencia.
```
#!/bin/bash
#bash nombre_shipt.sh

SEQS=/srv/home/chernandez/processing/bowtie
BIN=/srv/home/mromero/mbin/bowtie2

COUNT=0
for FAA in `ls *.fastq | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/bowtie2 -x calabacitagenome -1 $SEQS/$FAA"_R1.fastq" -2 $SEQS/$FAA"_R2.fastq" -S $FAA.sam"  -p 20 >>$*.$COUNT.scr
chmod +x *.scr; done
```
c) Se generan archivos en formato sam que tienen que convertirse a bam para continuar con el procesamiento
```
#!/bin/bash
#bash nombre_shipt.sh

FILES=/srv/home/chernandez/processing/bowtie
BIN=/srv/home/chernandez/binc/samtools-1.3.1

COUNT=0
for FAA in `ls *.sam | sed -e 's/.sam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/samtools view -bS $FILES/$FAA".sam"  > $FAA.bam"  >>$*.$COUNT.scr
chmod +x *.scr; done
```
d) Obtener las secuencias que no mapearon contra el hospedero.
```
#!/bin/bash
#bash nombre_shipt.sh

FILES=/srv/home/chernandez/processing/bowtie
BIN=/srv/home/chernandez/binc/samtools-1.3.1

COUNT=0
for FAA in `ls *.bam | sed -e 's/.bam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/samtools view -b -f 12 -F 256  $FILES/$FAA".bam"  > $FAA_"um.bam  >>$*.$COUNT.scr
chmod +x *.scr; done
```
e) Se orden las secuencias para que los pares (fordward y reverse) aparezcan juntas.
```
#!/bin/bash
#bash nombre_shipt.sh

FILES=/srv/home/chernandez/processing/bowtie
BIN=/srv/home/chernandez/binc/samtools-1.3.1

COUNT=0
for FAA in `ls *_um.bam | sed -e 's/_um.bam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/samtools sort -n  $FILES/$FAA"_um.bam"  > $FAA"_sorted.bam  >>$*.$COUNT.scr
chmod +x *.scr; done
```
f) Recuperar las secuencias forward y reverse por separado
```
#!/bin/bash
#bash nombre_shipt.sh

FILES=/srv/home/chernandez/processing/bowtie
BIN=/srv/home/chernandez/binc/bedtools2/bin

COUNT=0
for FAA in `ls *_sorted.bam  | sed -e 's/_sorted.bam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "$BIN/bedtools bamtofastq -i  $FILES/$FAA"_sorted.bam"  -fq $FILES/$FAA"_filtered_R1.fastq" -fq2 $FILES/$FAA"_filtered_R2.fastq >> $*.$COUNT.sh
chmod +x *.sh; done
>>>>>>> EdiciónCristóbal
```
1. Filtrado de calidad (Trimmomatic)
```
#!/bin/bash
#bash nombre_shipt.sh

SEQS=/srv/home/chernandez/processing/assembly
BIN=/srv/home/chernandez/binc/Trimmomatic-0.36

COUNT=0
for FAA in `ls *.fastq | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.sh
echo "java -jar $BIN/trimmomatic-0.36.jar PE -threads 2 -phred33 -trimlog trim.log $SEQS/$FAA"_filtered_R1.fastq" $SEQS/$FAA"_filtered_R2.fastq" $SEQS/$FAA"_paired_R1.fastq" $S
EQS/$FAA"_unpaired_R1.fastq" $SEQS/$FAA"_paired_R2.fastq" $SEQS/$FAA"_unpaired_R2.fastq" ILLUMINACLIP:$BIN/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:1
5 MINLEN:36" >>$*.$COUNT.scr
chmod +x *.scr; done
```
2. Ensamble (metaSPADES o megaHIT)
-Con metaSPADES:
```
#!/bin/bash

SEQS=/home/cristobal/spades
BIN=/home/cristobal/binc/SPAdes-3.12.0-Linux/bin

COUNT=0
for FAA in `ls *_paired_R*.fastq | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo  "$BIN"/spades.py --meta -1 "$SEQS/$FAA"_paired_R1.fastq -2 "$SEQS/$FAA"_paired_R2.fastq -s "$SEQS/$FAA"_all_unpaired.fastq -o "$SEQS/$FAA/contig_$FAA/ >>$*.$COUNT.scr
chmod +x *.scr; done
```
-Con megaHIT
```
#!/bin/bash
#bash nombre_shipt.sh <nombre-del-trabajo>

SEQS=/srv/home/chernandez/processing/assembly
OUT=/srv/home/chernandez/processing/assembly/megahit
BIN=/srv/home/mromero/mbin

COUNT=0
for FAA in `ls *_paired_R1.fastq | sed -e 's/_paired_R1.fastq//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "$BIN/megahit/megahit -m 24000000000 -t 40 -1 $SEQS/$FAA"_paired_R1.fastq" -2 $SEQ
S/$FAA"_paired_R2.fastq" -o $OUT/$FAA.megahit_result" >>$*.$COUNT.scr
```
3. Estadísticas de ensamblado. 
Se puede realizar para los contigs largos, cortos o el ensamble híbrido (contigs cortos + contigs largos):
```
python /home/miguel/mbin/quast-5.0.0/quast.py *contigs.fasta
```
4. Mapeo de _reads_ crudos _vs contigs_, se utiliza bbmap
```
#!/bin/bash

CONT=/home/cristobal/annotation/contigs
SEQS=/home/cristobal/annotation
BIN=/home/cristobal/binc/bbmap

COUNT=0
for FAA in `ls *paired_R1* | sed -e 's/_paired_R1.fastq//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "$BIN"/bbwrap.sh ref="$CONT/$FAA"_contigs.fa in="$SEQS/
$FAA"_paired_R1.fastq in2="$SEQS/$FAA"_paired_R2.fastq build="$COUNT" out="$SEQS/$FAA".sam kfil
ter=22 subfilter=15 maxindel=80 >> $*.$COUNT.scr
chmod +x *.scr; done
```
5. Obtener _reads_ no mapeados en el ensamblado
```
#!/bin/bash

SEQS=/home/cristobal/annotation

COUNT=0
for FAA in `ls *.sam | sed -e 's/.sam//g'`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo /home/cristobal/binc/samtools-1.2/samtools view -u -f 12 -F 256 "$SEQS/$FAA".sam "|" /home/cristobal/binc/samtools-1.2/samtools bam2fq - ">" "$FAA"_ump.fastq  >> $*.$COUNT
.scr
chmod +x *.scr; done
```
6. Segundo ensamble de _reads_ no mapeados con _Velvet_ con cobertura 2X
```
#!/bin/bash

SEQS=/home/cristobal/contigs
BIN=/home/cristobal/binc/velvet

COUNT=0
for FAA in `ls *_ump.fastq.gz | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo  "$BIN/velveth $SEQS/$FAA"_assemblyVelvet 21 -fastq -short "$SEQS/$FAA"_all_unpaired.fastq.gz -shortPaired "$SEQS/$FAA"_ump.fastq.gz >>$*.$COUNT.scr
echo "$BIN/velvetg $SEQS/$FAA"_assemblyVelvet -exp_cov 2 -ins_length 350 -min_contig_lgth 100 >> $*.$COUNT.scr
chmod +x *.scr; done
```
7. Mapeo de _reads_ crudos _vs_ ensamblado con _Velvet_
8. Descartar contigs de <100 pb; unir el ensamble del paso 2 contra el ensamble de 6
9. Predicción de ORFs con _Prodigal_
10. Anotación de los archivos de proteínas predichas de _Prodigal_ contra M5NR usando el API del RAST
11. Generar tabla de anotación de AfOTU (_annotated function OTU_).
12. Separar ORFs sin _hit_ al M5NR y generar un multifasta.
13. Clustering al 70% de las proteínas predichas sin hit al M5NR. Estos serán dominados como PUFO _(protein family of unknown origin)_
14. Fusionar la tabla del paso 11 con la  de los PUFOs del paso 13
15. Para la anotación, usar el SEED y los subsistemas para generar una tabla de taxonomía de descripción general; para la anotación fina usar las anotaciones de REFSEQ a nivel de funciones en otra tabla de taxonomía. Agregar a los PUFOs un identificador del estilo PUFO_1, PUFO_2, etc.
16. Cargar estos datos de phyloseq
