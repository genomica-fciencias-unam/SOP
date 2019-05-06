#!/bin/bash
#Genera un archivo con una lista de los integrantes de cada familia de genes/OTUs a partir de un archivo OTU
awk '{for(i=2;i<=NF;i++){printf "%s\n", $i}; printf "\n"}' fago.otu | perl -pe 'print "START\n" if $. == 1'| perl -pe 's/^\n/START\n/g' | awk '/START/{x='F'++i;next}{print >x".list";}'
