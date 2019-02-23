---
layout: default
title: Etiquetas para muestras e identificadores
---

Las etiquetas pensadas desde el inicio de los proyectos nos permiten tener trazabilidad al origen, tipo de muestra, colector y otros datos relevantes.

Para renombrar archivos de forma masiva se les recomienda utilizar el programa gprename, evaluar el nombre previo del nombre y proceder a renombrar. El comando gprename puede utilizarse de forma gráfica si se conectan a los servidores exportando las X o en sus sesiones de linux locales.

Los elementos mínimos de etiquetado de muestra son: 

1. Dos a tres iniciales en mayùsculas por persona (_i.e._ LA)
2. Clave del proyecto (tres letras)
	-PLA: Asociado a plantas/suelo
    -MET: Asociado a metro/humanos
    -ANI: Asociado a animales
    -AGU: Asociado a agua
3. Código identificador de muestra (tres dígitos: 000, 001, etc.).
4. Fecha de colecta DDMMAA
4. Còdigo identificador de tipo de procesamiento:
    -S: 16S
    -I: ITS
    -M: Metagenoma _shotgun_
    
Ejemplos:

LAPLA001-230219-S.R1.fastq.gz

 -Corresponde a la muestra 001 colectada por Luis Alcaraz, de un microbioma asociado a plantas, con datos del gen 16S rRNA. La muestra se obtuvo el 23 de febrero de 2019. 
 
 Los sufijos como fastq.gz son particulares de cada tipo de archivo, en este caso corresponde a la lectura 1 en formato comprimido fastq.
 
