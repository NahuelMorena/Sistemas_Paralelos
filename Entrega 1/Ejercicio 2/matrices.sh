#!/bin/bash

#Se compila el programa 
gcc -o matrices matrices.c -lm

#Se setean variables
array=(512 1024 2048 4096)

#Se crea el informe
> informe.txt
echo -e "Estudio de tiempos de ejecución del algoritmo de matrices con diferentes dimenciones\n\n" >> informe.txt


#Comienzan las pruebas
for size in ${array[@]}; do 
    let aux=$size/2
    echo "Prueba de ejecución con matrices de dimensión de ${size} y bloques de ${aux}"
    ./matrices ${size} ${aux} >> informe.txt
done

echo "Finalizan las pruebas, resultados en informe.txt"
