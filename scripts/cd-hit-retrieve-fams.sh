#!/bin/bash
#Genera un archivo con una lista de los integrantes de cada familia de genes/OTUs a partir de un archivo OTU (ver procedimientos 16S/ITS)
awk '{for(i=2;i<=NF;i++){printf "%s\n", $i}; printf "\n"}' fago.otu | perl -pe 'print "START\n" if $. == 1'| perl -pe 's/^\n/START\n/g' | awk '/START/{x='F'++i;next}{print >x".list";}'

##Extrae las secuencias y genera fastas individuales a partir del archivo secuencias.fasta
for N in `ls *.list`; do perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $N secuencias.fasta >$N.fasta; done
