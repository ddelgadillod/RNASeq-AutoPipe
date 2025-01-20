# RNASeq-AutoPipe
Flujo automático de procesamiento de datos de RNA Seq para obtener expresión diferencial de genes.

Este repositorio es una colección de Jupyter Notebooks y un preliminar de una alicación Web que pimplementan una automatización de un flujo de análisis de datos RNASeq, para obtener una lista de genes expresados diferencialmente.

## Prerequisitos
### Archivos entrada
- Lecturas RNASeq
- Genoma de referencia
- Tabla de conteos*
- Archivo de diseño o metadatos*
 *Suficientes para ejecutar la aplicación Web. 
### Software
1. Instalar [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) o [Mamba](https://github.com/conda-forge/miniforge), seguir las instrucciones seg[un sea de su preferencia.
2. Crear el ambiente a partir del archivo `rna-env.yml`:
```
conda env create -f rna-env.yml
```
3. Activar el ambiente:
```
conda activate RNASeq-Auto
```
5. Ejecutar Jupyter Notebook
6. Seguir la instrucciones y navegar para abrir los Notebook(Archivos de extensión `.ipynb).
## Notebooks
- RNASeq-ReadsQC-Aln.ipyn
- RNASeq-DataCheck.ipynb
- R-Analysis-PreDESeq2.ipynb
- DESeq2-Preprocessing.ipynb
- RNASeq-R-Analysis-PE-DESeq2.ipynb

## Aplicación Web
La aplicación Web ha sido desarrollada en `Shiny Python`, es una versión muy prematura con interfaz gráfica, que implementa el flujo propuesto por PyDESEq2 en sus tutoriales y permite obtener gráficas de forma automatica, a partir de la carga de la tabla de conteos y el archivo de metadatos o diseño. En el siguiente video se pude ver en forma resumida su funcionamiento.

[![Prueba de concepto interfaz gráfica del flujo automático](https://img.youtube.com/vi/Ic0vBIdgoXg/0.jpg)](https://www.youtube.com/watch?v=Ic0vBIdgoXg)


## Referencias
- BWA
- featureCounts
- Kallisto
- DESeq2
