# Mini curso abril 2021: Panorama general de análisis de datos de RNA-seq con R

Este es el material del mini curso "Panorama general de análisis de datos de RNA-seq con R" ofrecido por la Red Mexicana de Bioinformática.

Fecha: 16 de abril de 2021

Instructores:

- [M.C. Ana Beatriz Villaseñor Altamirano]()
- [M.C. Rodolfo Luis Chávez Domínguez]()

Para descargar este material da click en el botón **Code** y selecciona la opción **Download zip**. Si lo prefieres, puedes clonar el repositorio en tu computadora usando git clone desde tu terminal.

## Software

Para realizar la práctica del curso se requieren de los siguientes programas y paquetes:

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
* [Salmon](https://github.com/COMBINE-lab/salmon/releases) (opcional)
* R > 4.0
* RStudio
* Bioconductor >= 3.12, instalarlo con:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
```

* Paqueteria de R descargada con Bioconductor:
	* _tximeta_
	* _SummarizedExperiment_
	* _DESeq2_
	* _PCAtools_

Instalarlos con:

```
paquetes = c("DESeq2", "tximeta", "PCAtools", "SummarizedExperiment")
BiocManager::install(paquetes)
```

* Paqueteria de R descargada del repositorio CRAN
	* _tidyverse_
	* _pheatmap_

Instalarlos con:

```
install.packages()
```

## Archivos

Si decides realizar toda la práctica necesitarás estos archivos. Descargar el archivo en formato zip llamado **Files.zip** desde este [vínculo](https://drive.google.com/file/d/1n8dQHFotDc-fHxrLZOPVFoS8D29Pf5SZ/view?usp=sharing) y colocalo en el folder principal de este repositiorio. Descomprime el archivo y sigue las indicaciones de los instructores.
En el archivo **Files.zip** encontrarás:

* `Archivos_fastq/`: Folder que contiene los archivos en formato .fastq
* `Salmon_quants/`: Folder que contiene las matrices de cuentas generadas por Salmon
* `genecode.v37.transcripts.fa`: Transcriptoma de referencia para realizar el alineamiento y conteo con Salmon
* `ObjetosR.RData`: Archivo de R que contiene los objetos de *SummarizedExperiment* y *DESeq* para realizar el análisis de expresión diferencial. **Descarga este [objeto](https://drive.google.com/file/d/1_j2Py3ifDRN-O_6ygq627xuKTAI2HmGb/view?usp=sharing) si no tienes suficiente espacio en tu computadora** 

## Descripción general

|    | Tema                                                        | Práctica | Tiempo (min) | Encargado |
|----|-------------------------------------------------------------|--------|--------------|-----------|
| 1  | Introducción a RNA-seq                                      | NA     | 5            | AnaB      |
| 2  | Organización del proyecto (folder de datos, fig, scripts)   | NA     | 5            | Rodolfo   |
| 3  | Flujo de trabajo para el análisis de RNA-seq                | NA     | 5            | AnaB      |
| 4  | ¿Qué tipo de archivos se necesitan?                         | NA     | 5            | Rodolfo   |
| 5  | Limpieza de datos (Trimming y filtrado)                     | Si     | 10           | AnaB      |
| 6  | Alineadores (Salmon)                                        | Si     | 10           | Rodolfo   |
| 7  | tximeta, tximport                                           | Si     | 10           | AnaB      |
| 8  | Objeto SummarizedExperiment en R                            | Si     | 10           | Rodolfo   |
| 9  | Normalización de los datos (RPKM, TPM)                      | Si     | 15           | AnaB      |
| 10 | PCA (outliers)                                              | Si     | 15           | Rodolfo   |
| 11 | DE (analysis)                                               | Si     | 15           | AnaB      |
| 12 | DE (visualization/resultados)                               | Si     | 15           | Rodolfo   |



