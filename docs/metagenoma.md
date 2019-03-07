

Para los metagenomas:

Los pasos generales son: 

0. Mapeo de secuencias contra el hospedero (Bowtie2)
```
comando del mal
```
1. Filtrado de calidad (Trimmomatic)
2. Ensamble (metaSPADES o megaHIT)
3. Estadísticas de ensamblado
4. Mapeo de _reads_ crudos _vs contigs_ 
5. Obtener _reads_ no mapeados en el ensamblado
6. Segundo ensamble de _reads_ no mapeados con _Velvet_ con cobertura 2X
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
