
##Protocolo limpieza y procesado inicial de amplicones
   

#Visualizar la calidad de las secuencias y armar una tabla con la siguiente información:

       Cantidad de READS | Longitud promedio reads | Secuencias Ensambladas
Muestra


Usar FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) Almacenar los plots en html en el directorio QC/ 


```bash
ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq >lista #genera la lista
ln -s scritps/assembly.sh assembly.sh #link simbólico al script de ensamblado
bash assembly.sh nov1616assembly #se corre el script que genera los trabajos de ensamble usando PANDASEQ (requerimento previo)
for N in `ls *.scr`; do qsub $N; done # se mandan al cluster los trabajos de ensamblado
```

El script de pandaseq utiliza los siguientes parámetros:

Particularly the following options were modified:

−o minoverlap
Sets the minimum overlap between forward and reverse reads. By default, this is at least one nucleotide of overlap. Raising this number does not generally increase the quality of the output as alignments with small overlaps tend to score poorly and are discarded anyway. This was chosen to be a value of 15, rather than 1. 
−O maxoverlap
Sets the maximum overlap between forward and reverse reads. By default, this is the read length. In highly overlapping sequences (i.e., those where the end of one read precede the start of the other), this parameter should be set to the sum of the input reads, or a value larger than that. The maxoverlap was set to default which is the read length
−t threshold
The score, between zero and one, that a sequence must meet to be kept in the output. Any alignments lower than this will be discarded as low quality. Increasing this number will not necessarily prevent uncalled bases (Ns) from appearing in the final sequence. It is also used as the threshold to match primers, if primers are supplied. The default value is 0.6. Our chosen value was 0.95 in a 0-1 scale. 

Length was filtered up to a size of 470 bp because the expected amplicon size is around 460 bp from the 16S rRNA gene V3-V4 primers (http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf).


```bash
#renombrar las secuencias con identificadores con el nombre de la secuencia y numerarlas:
#>Secuencia_0 


perl scripts/header.fasta.numbers.pl PREFIX nombre_del_archivo.fasta 
#Hacer este paso por cada muestra. Esto sirve para tener etiquetas individuales por muestra
```