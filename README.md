drugAnnotation
==============

Annotate drugs based on SMILE string or drug names from ChemSpider, ChEMBL, ChEBI, KEGG, clinicaltrials.gov.

Create ontology with functional terms from drug annotation and perform tests of functional enrichment of drugs with a Fisher exact test.

# Usage
`python drugAnnotation.py [OPTIONS] <file.csv> <output directory>`

and

`python drugEnrichment.py [OPTIONS] test.csv universe.csv <output directory>`



# Requirements
To install all requirements simply run:
`pip install -r requirements.txt`

# Drug annotation
The main goal is the functional annotation of drugs from different sources, so I extract only information that allows easier identification of the drugs (IDs, SMILES, etc...) and any information related with the drug function.

Annotation from ChemSpider, ChEMBL and ChEBI is performed based on SMILE structure.
KEGG and clinicaltrials.gov can currently (to my knowledge) only be queried by drug name, which is obviously prone to higher error.

By default annotation is only performed to services that allow querying by SMILE.

## ChemSpider
Match is made based on SMILES.

Set up you access token with the environment variable `CHEMSPIDER_SECURITY_TOKEN` in bash like this: `export $CHEMSPIDER_SECURITY_TOKEN=<your chemspider token>`. You can also add it to your `.bashrc` or `.bash_profile` for convinience.

## ChEMBL
Match is made based on SMILES.

## ChEBI
Match is made based on SMILES with at least 70% identity.

## KEGG
Name-based. Only KEGG ids and names are annotated.

## clinicaltrials.gov
Only completed studies are inspected. Study title, url and outcome are extracted.

## Functional terms ontology
Based on ChEBI terms, an ontology is built for all drugs with annotation, which allows testing for term enrichment in drug subsets.


## Output
Individual csv files are created for each origin of annotation.
A concatenated set is in `drugAnnotation.csv`.
All files produced are in `/results`.


# Test the functional enrichment of subsets of drugs

Once with drugs annotated with terms from the chebi ontology, it is possible to test for over-representation of ontology terms with a Fisher's exact test.

Run the enrichment test part with the following command: `python drugEnrichment.py [OPTIONS] test.csv universe.csv`, where `test.csv` is a csv file with a subset of drugs previously annotated with `drugAnnotation.py` and `universe.csv` is all of the drugs previously annotated.

By default, test results for all ontology terms tested (significant or not) will be outputed, leaving to the user the liberty of setting a threashold based on a p-value (this way, even unsignificant terms can be explored).

## Using a score for each drug used in chemical screens

When running chemical screens, the effect of interesting drugs may have contrasting biological outcomes (e.g. improved cell death or survival). To test grops of drugs with such oposing characteristics separately, just indicate a column in the csv file with a score for this metric with the `-s` or `--score-column` argument. Positive and negative values will be analysed separately.


# Running tests

`python drugAnnotation.py testAnnotation.csv`

`python drugEnrichment.py --id-column-test 3 testEnrichment.csv annotation/universeOntology.csv`
