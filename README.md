# Bioinformática: Trabajo práctico


## Intro
Genebank source: https://www.ncbi.nlm.nih.gov/nuccore/NG_012386.1?from=5001&to=58286&report=genbank

## Requerimientos

### Python 3.6
Para ejecutar el proyecto es necesario tener instalado la versión 3.6 de Python.

Ejecutar el archivo `init.sh` para instalar el entorno virtual y sus dependencias.

## Ejercicios

Para ejecutar cualquier ejercicio es tan simple como correr:
```
./venv/bin/python src/exX.py
```

Todos los ejercicios tienen el comando `-h` que imprime la ayuda de el mismo.

Se incluyeron una serie de scripts que tienen invocaciones por defecto utilizando los datos provistos.
Para su correcto funcionamiento deben de ser ejecutados en orden (ejercicio 1 primero, el 2 segundo y así sucesivamente).

### Ejercicio 1
Para ejecutar el ejercicio 1:
```
./ej1.sh
```
El script convertirá un archivo en formato GeneBank a format Fasta para ser utilizado en el resto de los ejercicios.

### Ejercicio 2

Si se quiere ejecutar el ejercicio 2 utilizando una base de datos BLAST local, se deberá contar con:

- El comando `blastx` en el `PATH` del sistema operativo
- La base de datos `swissprot` correctamente instalada y debe de poder ser encontrada por el comando `blastx`

Para configurar correctamente BLAST, referirse a la [documentación oficial](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

Es posible ejecutar el ejercicio 2 utilizando un BLAST local o remoto (ver ayuda para instrucciones sobre su configuración)

Para ejecutar una versión preestablecida utilizando BLAST remoto ejecutar el siguiente script:
```
./ej2.sh
```
