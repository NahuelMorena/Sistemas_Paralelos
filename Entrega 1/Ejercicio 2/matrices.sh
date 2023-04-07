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
    echo "Prueba de ejecución con matrices de dimeción de ${size}"
    ./matrices ${size} >> informe.txt
done

echo "Finalizan las pruebas, resultados en informe.txt"
