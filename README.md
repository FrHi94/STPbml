# STPbml (Spatial Transcriptomics of *Plasmodium berghei* infected mouse liver)

This repository contains code and data required to reproduce data described in the study "Host-Pathogen Interactions in the Plasmodium-Infected Mouse Liver at Spatial and Single-Cell Resolution" currently available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.12.22.573046v1). Larger datasets are available in a seperate zenodo repository: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8386528.svg)](https://doi.org/10.5281/zenodo.8386528). 

## system requirements 

All code can be reproduced using R and python software. We recommend using R ```version 4.0.5``` and python ```3.10.8``` or higher. Moreover to reprodcue the generation of distance measurements the hepaquery package is required. Installation instructionr, requirements and dependencies for this package can be found [here](https://github.com/almaan/ST-mLiver). 

## content

- ```data```

  - ```annotations``` 

    - contains GFF (General Feature Format) file for aligment of sequencing reads to joined *Plasmodium berghei* and *Mus musculus* genome using [stpipeline](https://github.com/SpatialTranscriptomicsResearch/st_pipeline) : *Mus_musculus.GRCm38.101_PlasmoDB-48_PbergheiANKA.gff.gz*

    - gene annotations including corresponding Gene Ontology terms for joined Genome: *20230222_annotation_Symbol_GO.txt*

    - gene annoations including gene symbols, biotype etc. for *Mus musculus* : *mus_pb_annotation_seurat.tsv*

    - gene annoations including gene symbols, biotype etc. for *Plasmodium berghei*: *pb_annotation_seurat.tsv*

  - ```celltype_proportions``` contains an output folder for each ST and Visium section with cell type proportion results of the [stereoscope](https://github.com/almaan/stereoscope) cell type deconvolotion and corresponding log files for the analysis. 

  - ```cluster``` contains cluster information for each spot for each individual sample of the ST analysis 

  - ```distances``` contains distance calculations performed using the hepaquery package and the yaml files in the ```yaml``` files folder and generated masks available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8386528.svg)](https://doi.org/10.5281/zenodo.8386528). 

  - ```gene_shortlists``` contains curated gene lists (by literature research) from expression-by-distance correlation analysis for parasite distances ordered by timepoint. 


- ```res```

  - contains folder to download seurat objects from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8386528.svg)](https://doi.org/10.5281/zenodo.8386528) to rerun analysis and generate respective plots and data shown in the study. 


- ```scripts```

 - ```R``` contains R code scripts (Notebooks in .RMD format or regulary .R scripts)

 - ```python``` contains python code scripts (Notebooks in .ipynb format)
