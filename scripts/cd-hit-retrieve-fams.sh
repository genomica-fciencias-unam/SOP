#!/bin/bash
#Luis D. Alcaraz mayo 2019

#Genera un archivo con una lista de los integrantes de cada familia de genes/OTUs a partir de un archivo OTU (ver procedimientos 16S/ITS)
awk '{for(i=2;i<=NF;i++){printf "%s\n", $i}; printf "\n"}' fago.otu | perl -pe 'print "START\n" if $. == 1'| perl -pe 's/^\n/START\n/g' | awk '/START/{x='F'++i;next}{print >x".list";}'

##Extrae las secuencias y genera fastas individuales a partir del archivo secuencias.fasta
for N in `ls *.list`; do perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $N secuencias.fasta >$N.fasta; done

#Cuenta el número de secuencias por familia y excluye a las que son de 1 secuencia (singletons)
grep -c ">" *.fasta | egrep ":1$" | perl -pe 's/\:1$//g' >single
grep -v -f single all  > align #genera la lista de archivos a alinear

#Genera un link simbólico con las secuencias singletons para que tengan un archivo .aln (de "alineamiento"

for L in `cat single`; do echo "ln -s $L $L.aln"; done | bash

#Alinea las secuencias (en máquina sencilla):
#for A in `cat align`; do echo "muscle -in $A -out $A.aln"; done | bash

#EN deep-thought / cuallicua (manda 1 alineamiento por procesador)
for A in `cat align`; do qsub -N NOMBRE_TRABAJO -b y -j y -cwd -V "muscle -in $A -out $A.aln"; done 


#Genera hmms con cada archivo de alineamiento

for N in `ls *.aln`; do hmmbuild $N.hmm $N; done


#EOF

