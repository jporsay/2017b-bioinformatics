# Bioinformática: Trabajo práctico


## Introducción

En el siguiente trabajo realizaremos análisis sobre la Tuberculosis Sclerosa, específicamente sobre el gen **TSC1**, gen que se encuentra en el cromosoma 9 y es encargado de encodear la proteína llamada Hamartina.

Fuente para el archivo genebank: [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NG_012386.1?from=5001&to=58286&report=genbank)

## Requerimientos

### Python
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
El script convertirá un archivo en formato GeneBank a formato Fasta para ser utilizado en el resto de los ejercicios.

### Ejercicio 2
#### Ejercicio 2.a

Si se quiere ejecutar el ejercicio 2 utilizando una base de datos BLAST local, se deberá contar con:

- El comando `blastx` en el `PATH` del sistema operativo
- La base de datos `swissprot` correctamente instalada y debe de poder ser encontrada por el comando `blastx`

Para configurar correctamente BLAST, referirse a la [documentación oficial](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

Es posible ejecutar el ejercicio 2 utilizando un BLAST local o remoto (ver ayuda para instrucciones sobre su configuración)

Para ejecutar una versión preestablecida utilizando BLAST remoto ejecutar el siguiente script:
```
./ej2.sh
```

#### Ejercicio 2.b

Podemos ver que los primeros 10 alineamientos están dominados por matches con otras proteínas humanas, sin embargo podemos también ver que aparecen los genes de TSC1 de ratas y ratones con un *Expect value* **muy** cercano a cero, lo que significa que el match fue muy significativo.

A su vez, podemos ver que uno de los matches pertenece a un organismo denominado *PANTR*, o *Pan troglodytes* (Chimpanzee).

Ambos resultados son interesantes ya que muestran las similitudes en ciertas proteínas que compartimos con los chimpanzees, ratas y ratones.


### Ejercicio 3

#### Multiple Sequence Alignment
El ejercicio 3 fue resuelto utilizando un servicio externo para realizar el MSA: https://www.ebi.ac.uk/Tools/msa/mafft/

Analizando el archivo generado en el ejercicio 2, *ej2_out.txt*, obtuvimos los **id** de los 10 mejores alineamientos y, utilizando estos IDs, realizamos consultas a la base de datos de BLAST para obtener las secuencias mediante el siguiente comando:

```
blastdbcmd -db swissprot -entry ID
```

Con cada salida de ese programa, generamos un archivo llamado **msa.fasta**, el cual fue utilizado en el servicio mencionado arriba para realizar el alineamiento múltiple.
La salida se encuentra en **aln-fasta.fasta**.

El script `ex3.py` es un punto de entrada para trabajar con la librería `AlignIO` de `BioPython` para manipular el MSA. Actualmente consume el archivo `msa.fasta`.
