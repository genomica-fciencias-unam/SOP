#!/bin/bash
#uso:
#./cd-clust.sh archivo.fas #corre hasta 100 procesos paralelos que se pueden modificar con los calificadores --S (número de trabajos en los que se parte) y --Q (número de trabajos empleados), la segunda parte reinicia la cola en cada iteración hasta que finaliza el trabajo, se puede modificar el programa (cd-hit, cd-hit-est) 


cd-hit-para.pl -i $* -o $*out --P cd-hit-est -c 0.97 -s 0.97 --S 100 --Q 100 --T "SGE" 

until cd-hit-para.pl -i $* -o $*out --P cd-hit-est -c 0.97 -s 0.97 --S 100 --Q 100 --T "SGE" --R $*out.restart | grep -q "CPU"; do echo "siguiente iteración"; done

