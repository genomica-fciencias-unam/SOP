# Correr promer en un ambiente SGE
El script promer de mummer4 tiene una función que automáticamente detecta el número de núcleos disponibles para paralelizar la ejecución de mgaps, lo cual no siempre es deseable.
Para definir el número de núcleos a usar hay que usar la variable "OMP_NUM_THREADS" en la terminal en donde se va a correr promer y asignar el mismo número de núcleos al ejecutar el qsub.

Ejemplo de script a correr: 

script "promer.scr"
```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
set OMP_NUM_THREADS=30
export OMP_NUM_THREADS=30
/home/mromero/mbin/mummer4/bin/promer -p prefix reference.fas query.fas
```

Ejemplo de llamado al SGE:

```
qsub -pe completenode 30 -l h_vmem=150G promer.scr
```
