drugAnnotation
==============

Annotate drugs based on SMILE string or drug names from ChemSpider, ChEMBL, ChEBI, KEGG, clinicaltrials.gov and create ontology with functional terms from drug annotation.

# Usage
usage: `python drugAnnotation.py [OPTIONS] file.csv`

positional arguments:
  infile                CSV file with drug smiles (and optionally names).

optional arguments:
  -h, --help            show this help message and exit
  -s SMILESCOLUMN, --smiles-column SMILESCOLUMN
                        Specify the column with the smile strings. Default is
                        column 1 (1-based).
  -n, --query-name      Annotate drugs by name as well as by smile. Off by
                        default.
  --names-column NAMESCOLUMN
                        Specify the column with drug names if available.
                        Default is column 2 (1-based).
  -l LOGFILE, --logfile LOGFILE
                        Specify the name of the log file.
  -v, --verbose         Verbose behaviour. Will print information to stdout.


# Drug annotation
The main goal is the functional annotation of drugs from different sources, so I extract only information that allows easier identification of the drugs (IDs, SMILES, etc...) and any information related with the drug function.

Annotation from ChemSpider, ChEMBL and ChEBI is performed based on SMILE structure.
KEGG and clinicaltrials.gov can currently (to my knowledge) only be queried by drug name, which is obviously prone to higher error.

By default annotation is only performed to services that allow querying by SMILE.

## ChemSpider
Match is made based on SMILES.

## ChEMBL
Match is made based on SMILES.

## ChEBI
Match is made based on SMILES with at least 70% identity.

## KEGG
Name-based. Only KEGG ids and names are annotated.

## clinicaltrials.gov
Only completed studies are inspected. Study title, url and outcome are extracted.

# Functional terms ontology
Based on ChEBI terms, an ontology is built for all drugs with annotation, which allows testing for term enrichment in drug subsets.
