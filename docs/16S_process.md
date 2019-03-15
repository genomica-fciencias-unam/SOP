
## Protocolo limpieza y unión de lecturas pareadas para amplicones
   

### Visualizar la calidad de las secuencias y armar una tabla con la siguiente información:

Usar [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) Almacenar los plots en html en el directorio QC/ 

Analizar la calidad de las secuencias y calcular las estadísticas de secuenciación crudas, generalmente viene en los reportes de Cuernavaca, Irapuato y Korea, aquí pueden reconfirmar. 


## Ensamble por CASPER **usar por defecto


Guardar los datos del % de lecturas ensambladas, ejemplo:


  -  Total number of reads     :      78186
  -  Number of merged reads    :      58986 (75.44%)
  -  Number of unmerged reads  :      19200 (24.56%)
  -  TIME for total processing :      3.247 sec

Armar un archivo tabular con la siguiente información:

Muestra | Cantidad de READS | Longitud promedio reads | # Secuencias Pareadas Unidas | Longitud promedio Seqs PE Unidas



```bash
ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq >lista #genera la lista
ln -s scritps/assemblyCASPER.sh assemblyCASPER.sh #link simbólico al script de ensamblado por PANDASEQ
bash assemblyPANDA.sh NOMBRE_TRABAJO #se corre el script que genera los trabajos de ensamble usando PANDASEQ (requerimento previo)
for N in `ls *.scr`; do qsub $N; done # se mandan al cluster los trabajos de ensamblado
```

## También se puede experimentar con PANDASEQ, hay que hacer cortado manual de las bases dependiendo de la calidad, este protocolo fue útil con MiSEQ con un 2x250 

```bash
#para ensamblar con PANDASEQ
ls *.fastq.gz | perl -pe 's/_R.*.fastq.gz//g' | sort | uniq >lista #genera la lista
ln -s scritps/assemblyPANDA.sh assemblyPANDA.sh #link simbólico al script de ensamblado por PANDASEQ
bash assemblyPANDA.sh NOMBRE_TRABAJO #se corre el script que genera los trabajos de ensamble usando PANDASEQ (requerimento previo)
for N in `ls *.scr`; do qsub $N; done # se mandan al cluster los trabajos de ensamblado
```

## Renombrar las secuencias a identificadores sencillos (ver sección de etiquetas)
renombrar las secuencias con identificadores con el nombre de la secuencia y numerarlas:
```
>Secuencia_0 
ATTAC
>Secuencia_1
TTCAT
>Secuencia_2
CCTAT
```

```bash
#Hacer este paso por cada muestra. Esto sirve para tener etiquetas individuales por muestra
perl scripts/header.fasta.numbers.pl PREFIX nombre_del_archivo.fasta 

#Concatenar todas las muestras del estudio

cat *.numbered.fasta >estudio_completo.fas

```

## Agrupamiento de secuencias al 97% de identidad

Si estás en el cluster _deep thought_ hay que utilizar el script de paralelización

Si estás usando el cluster _cuallicua_ se utiliza directo el script, pero se envia por el gestor de colas (qsub)


```bash
#En deep-thought:
ln -s scripts/cd-clust.sh . 
bash cd-clust.sh estudio_completo.fas

#En cuallicua:
qsub -N NOMBRE_TRABAJO -b y -j y -cwd -V "cd-hit-est -c 0.97 -T 25 -M 0 -i estudio_completo.fas -o output.clstr"
```


```bash
#Generar tabla de OTUs

perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g; s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' output.clstr.clstr >estudio.otus

#Obtener secuencias representativas
pick_rep_set.py -i estudio.otu -f estudio_completo.fas -o rep_set.fna

```

## Filtrado de secuencias que no son 16S (se puede adaptar también para ITS)

Primero hacer un blast contra una DB pequeña de 16S (OTUs 70% id)
Identificar los OTUs que dieron aciertos
Generar un archivo de trabajo libre de secuencias que no sean 16S


```bash
parallel_assign_taxonomy_blast.py -i rep_set1.fna -o no16S_screen -r /qiime/gg_otus-13_8-release/rep_set/70_otus.fasta -t /qiime/gg_otus-13_8-release/taxonomy/70_otu_taxonomy.txt

cat no16S_screen/repset_tax_assignments.txt | grep -c "No blast hit"

cat no16S_screen/repset_tax_assignments.txt | grep -v "No blast hit" | cut -f1 >ids_screened.txt

cat no16S_screen/repset_tax_assignments.txt | grep "No blast hit" | cut -f1 >ids_REMOVE_biom.txt


perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids_screened.txt rep_set.fna >rep_set.screened.fna #extrae las secuencias con match a 16S y hace un nuevo archivo representativo


#asignación taxonómica del 16S del dataset completo al 97%
parallel_assign_taxonomy_blast.py -i rep_set.screened.fna -o taxonomy -r /qiime/gg_otus-13_8-release/rep_set/97_otus.fasta -t /qiime/qiime_software/gg_otus-13_8-release/taxonomy/97_otu_taxonomy.txt





```

## Armar tabla de OTUs


```bash
make_otu_table.py -i estudio.otu -t taxonomy/taxonomy_repset_tax_assignments.txt -o estudio.biom 

#Quitar los OTUs sin hits de 16S y singletons
filter_otus_from_otu_table.py -i estudio.biom -e ids_REMOVE_biom.txt -o estudio_screened.biom -n2 ; mv estudio_screened.biom estudio.biom



#alinear secuencias para identificar químeras:
parallel_align_seqs_pynast.py -i rep_set.screened.fna -o chimalign -X estudio

#identificar quimeras
parallel_identify_chimeric_seqs.py -i chimalign/estudio.repset_aligned.fasta -r /qiime/gg_otus-13_8-release/rep_set_aligned/85_otus.fasta -o estudio.chimera.txt
awk '{print $1}' estudio.chimera.txt >ids_CHIMERA.txt


#eliminar quimeras del archivo biom

filter_otus_from_otu_table.py -i estudio.biom -e ids_CHIMERA.txt -o estudio_chimera.biom; mv estudio_chimera.biom estudio.biom



#generar biom tabular
biom convert --to-tsv -i estudio.biom -o estudio.biom.tsv --table-type "Taxon table" --header-key=taxonomy


```
