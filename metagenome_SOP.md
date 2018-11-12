
#Standard operative procedures for WGS metagenomics @ LDA lab
lalcaraz@ciencias.unam.mx


##Assembly and QC


```bash
#!/bin/bash
#guardar este script como quieras nombrarlo y ejecutarlo como se indica a continuación:
#bash nombre_script.sh <nombre-del-trabajo>
#ejecutarlo en un directorio con los archivos crudos de sequenciación (SEQ)
#LDA feb 2017

SEQS=/home/luis/WGS_metagenomas2017/raw_fastq
SALIDAS=/home/luis/WGS_metagenomas2017/raw_fastq
BIN=/home/luis/bin/Trimmomatic-0.36
BIN2=/home/luis/bin/velvet
BIN3=/home/luis/bin/MetaVelvet

COUNT=0
for FAA in `ls *.gz | perl -pe 's/\_.*//g' | sort | uniq`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "#$ -l h_vmem=12G" >>$*.$COUNT.scr


#Trimmomatic
echo "java -jar $BIN/trimmomatic-0.36.jar PE -threads 2 -phred33 -trimlog trim.log $SEQS/$FAA"_1.fastq.gz" $SEQS/$FAA"_2.fastq.gz" $SEQS/$FAA"paired_1.fastq.gz" $SEQS/$FAA"unpaired_1.fastq.gz" $SEQS/$FAA"paired_2.fastq.gz" $SEQS/$FAA"unpaired_2.fastq.gz" ILLUMINACLIP:$BIN/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" >>$*.$COUNT.scr


#velveth y velvetg
echo  "$BIN2/velveth $SEQS/$FAA"_assemblyV3" 21 -shortPaired -fastq $SEQS/$FAA"paired_1.fastq.gz" $SEQS/$FAA"paired_2.fastq.gz"" >>$*.$COUNT.scr 
echo "$BIN2/velvetg $SEQS/$FAA"_assemblyV3" -exp_cov 2 -ins_length 350" >>$*.$COUNT.scr

#metavelvet
#echo "$BIN3/meta-velvetg $SEQS/$FAA"_assembly" -ins_length 350"  >>$*.$COUNT.scr

done
```


```bash
qsub nombre-del-trabajo.scr
```


```bash
#!/usr/bin/perl
#remove small contigs; por ejemplo de al menos 100 pb:
#for N in `ls */*.fa`; do echo "perl small.pl 100 $N > $N"l100.fa"" ; done | bash
use strict;
use warnings;
 
my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}

```

##Gene prediction & annotation


```bash
#Prodigal descarga y parámetros para correr
wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux
chmod +x prodigal.linux 
for N in `ls *.fas`; do echo "qsub -N $N.scr -b y -j y -cwd -V \"/home/luis/bin/prodigal.linux -a $N.faa -p meta -i $N\""; done | bash

```


```bash
#!/bin/bash
#Se parte el fasta resultante entre el número de trabajos que se quiere enviar
#prc es el número de trabajos que se quiere enviar, en este cluster tenemos 40 núcleos,
#pero si le pedimos que use dos procesadores por trabajo lo dividimos en 20, modificar el valor
#prc de acuerdo a necesidades
#ejecutar: $~/bash nombre_script | bash  

seqs=`grep -c ">" $*`
prc=20
parte=$((seqs/prc))
echo "perl /usr/local/bin/partefasta $parte $*"


```


```bash
#En este bucle se mandan los trabajos a la cola para correr el programa eggnog-mapper
#utilizando la base de datos de bacterias, se piden 2 procesadores por proceso y se mandan al 
#cluster. Puedes cambiar el *.fas por cualquier otro wildcard. Es muy importante corregir
#los directorios de salida a un directorio donde tengan permisos de escritura. 

for N in `ls *.fas`; do echo "qsub -pe completenode 2 -N $N -b y -j y -cwd -V "python /databases/eggnog-mapper/emapper.py -i /databases/test/$N -o /databases/test/out$N -d bact --cpu 2""; done  | bash

#Cuando terminan los trabajos hay que concatenar las salidas:

cat *.annotations | grep -v "#" >outcca.annotation
cat *.hmm_hits | grep -v "#" >outcca.hmm_hit
cat *.seed_orthologs | grep -v "#" >outcca.seed_ortholog

#La explicación de los reultados la pueden encontrar aquí: https://github.com/jhcepas/eggnog-mapper/wiki/Results-Interpretation#project_nameemapperhmm_hits-file

```


```bash
#anotación por m5nr (local) en deep-though
#las bases de datos se encuentran disponibles en deep-thought en /databases

#ejermplo: blastp -db /databases/m5nr16may17/m5nr16may17 -word_size 6 -query /databases/test/cca.faa.0.fas -outfmt 6 -evalue 1e-10 -num_alignments 10 -num_threads 2 -out /databases/test/blastcca.faa.0.fas.bout
#para correr el comando en todos los nodos de deep-thought, con archivos de entrada terminados en *.fas:

#Correr con diamond!!

for N in `ls cca*.fas`; do echo "qsub -pe completenode 2 -N $N -b y -j y -cwd -V \"/home/luis/diamond blastp -p 2 -d /databases/m5nr16may17/m5nr.dmnd -q $N --outfmt 6 -e 1e-10 -o diamond$N.bout"


#for N in `ls cca*.fas`; do echo "qsub -pe completenode 2 -N $N -b y -j y -cwd -V \"blastp -db /databases/m5nr16may17/m5nr16may17 -word_size 6 -query /databases/test/$N -outfmt 6 -evalue 1e-10 -num_alignments 10 -num_threads 2 -out /databases/test/blast$N.bout\""; done | bash
#ten cuidado de revisar las rutas completas de ejecución, bases de datos y directorios de salida. 

#oneliner para concatenar salidas de blast, ordenarlas por valor de bitscore, remover duplicados con el mismo valor de bitscore se guardan en el archivo: best_uniq
cat *.bout | perl -pe ' $name_col=0; $score_col=11; while(<>) { s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) { push @names, $n }; if (! exists($max{$n}) || $s > $max{$n}) { $max{$n} = $s; $best{$n} = () }; if ($s == $max{$n}) { $best{$n} .= "$_\n" }; } for $n (@names) { print $best{$n} } ' >best;  perl -e ' $column=0; $unique=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' best >best_uniq; rm best
awk '{print $1"\t"$3"\t"$11"\t"$12"\t"$2}' best_uniq.tsv >best.simple.tsv #ligero reacomodo de columnas para la anotación en los siguientes pasos. 

cut -f2 best_uniq >test2 #identificadores de md5 para buscar en RAST-
/databases/SUBSYSTEMS.tsv #este archivo se obtuvo del api del RAST descargando todos los subsistemas, se filtro y edito manualmente
# http://api.metagenomics.anl.gov/1/m5nr/ontology?source=Subsystems&min_level=function obtener todo el json de la 
1285  perl -e ' $col1=9; $col2=0; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' simple+genbank+subsystems SUBSYSTEMS.tsv >complete.tsv
 
#editar la salida y quitar redundancia de md5 y otros identificadores se sugiere el siguiente orden:
#protein	%id	E-value	bitscore	md5	Accesion	Genbank annotation	Best_hit	Subsystem	Annotation subsystems	Level1	Level2	Level3	Level4
#cca_assembly_6	100	1E-113	336	906aad964cdabaf371a78677aef566eb	EAY63195.1	Primosomal protein n	Burkholderia cenocepacia PC184	SS23169	Helicase PriA essential for oriC/DnaA-independent DNA replication	DNA Metabolism	DNA replication	DNA-replication	Helicase PriA essential for oriC/DnaA-independent DNA replication

#sources M5NR http://api.metagenomics.anl.gov/m5nr/sources








```


```bash
#!/bin/bash
#para correr múltiples anotaciones en varios archivos hice este script (anota.sh) :
for N in `ls *.bout`; 
do
/databases/myM5NR/bin/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source RefSeq --md5 $N.best_uniq.busca | perl -pe ' $column=1; $unique=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' >$N.annotation_refseq;
/databases/myM5NR/bin/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source Subsystems --md5 $N.best_uniq.busca| perl -pe ' $column=1; $unique=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' >$N.annotation_Subsystems;
/databases/myM5NR/bin/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source SEED --md5 $N.best_uniq.busca | perl -pe ' $column=1; $unique=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' >$N.annotation_Seed;
perl -e ' $col1=1; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' $N.best_uniq $N.annotation_refseq  | cut -f1,2,3,10,11,13,15,16 >$N.hits+seed
perl -e ' $col1=1; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' $N.hits+seed $N.annotation_Subsystems| cut -f1-9 >$N.seed+subsystems;
#revisar
perl -e ' $col1=8; $col2=0; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' $N.seed+subsystems SUBSYSTEMS.tsv | cut -f1-8,10-14  >$N.seed+refseq+subsystems.tsv; 

done
```


```bash
#!/bin/bash
##spliter para archivos multifastas en formato de >secuencia_1 >seq_3 toma los nombres de las muestras y separa a multifastas individuales (pasar de repset o concatenados a archivos individuales)
echo obtaining sequence tags of $*
grep '>' $* | perl -pe 's/\_.*$//g;s/\>//g' | uniq >$*.samples

echo "/^$/    { next; }" >splitter.awk
echo "/^>/    { filename=\"$*.WARNING.notmatched\";}" >>splitter.awk

for SUBMIT in `cat $*.samples`
do

echo "/^>$SUBMIT/   { filename=\"$SUBMIT.fna\"; }" >>splitter.awk
    

done

echo "{ print >> filename;}" >>splitter.awk

# cat splitter.awk

awk -f splitter.awk $*

echo "this is all folks"

#EOF

```


```bash
#anotaciòn por m5nr patung 
#hay que generar los archivos condor y los scripts para hacerlos correr en el sistema de colas condor. 
#usar un archivo header que contenga solamente: #!/bin/bash

for N in `ls *.fas`; do echo "executable = $N.sh"  >$N.condor; echo "output = $N.out" >>$N.condor; echo "error =$N.out"  >>$N.condor; echo "log = $N.out" >>$N.condor; echo "request_cpus = 2" >>$N.condor; echo "queue 1" >>$N.condor; chmod +x *.condor; done

for N in `ls *.fas`; do cat header >$N.sh  ; echo "date" >>$N.sh; echo "/srv/home/lalcaraz/bin/blastp -query /srv/home/lalcaraz/shama/$N -db /srv/home/lalcaraz/m5nr16may17/m5nr16may17 -word_size 6 -outfmt 6 -evalue 1e-10 -num_alignments 10 -num_threads 2 -out /srv/home/lalcaraz/shama/$N.bout" >>$N.sh ; echo "date" >>$N.sh; chmod +x *.sh ; done

```

##Taxonomic binning & assignment


```bash
#kraken estas instrucciones son para correrlo en el cluster de Lancis patung

for N in `ls *.fas`; do echo "executable = $N.sh"  >$N.condor; echo "output = $N.out" >>$N.condor; echo "error =$N.out"  >>$N.condor; echo "log = $N.out" >>$N.condor; echo "request_cpus = 20" >>$N.condor; echo "queue 1" >>$N.condor; chmod +x *.condor; done

for N in `ls *.fas`; do cat header >$N.sh  ; echo "date" >>$N.sh; echo "/srv/home/lalcaraz/bin/kraken --db /srv/home/lalcaraz/ --fasta_input $N --threads 40 --output $N.kraken" >>$N.sh ;echo "date">>$N.sh; chmod +x *.sh ; done

#El archivo header solo contiene un: #!/bin/bash

for N in `ls *.kraken`; do kraken-translate --mpa-format --db /srv/home/lalcaraz/ $N >tax.$N; done

#todo en uno
for N in `ls tax.*`;do cat $N | perl -pe 's/_cov_/\t/g;s/\|/\t/g' >$N.tmp;  awk 'FNR==NR{if(m<NF)m=NF;next}{for(i=NF;i<m;i++)$(i+1)="\tNA"}1' $N.tmp $N.tmp >$N.parsed; rm $N.tmp; done


#por partes, si ya hiciste el último paso es la unión de los dos siguientes:
cat tax.$N.kraken  | perl -pe's/_cov_/\t/g;s/\|/\t/g' >test #en este caso era un ensamblado de velvet, se utiliza la cobertura como una métrica y se pone como una columna separada
awk 'FNR==NR{if(m<NF)m=NF;next}{for(i=NF;i<m;i++)$(i+1)="\tNA"}1' test test   >split.txt #en esta linea se evalua la cantidad de columnas máximas, la salida no es homogenea, si no encuentra el número igual rellena con NAs separados por tabuladores


```
