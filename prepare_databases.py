#!/usr/bin/env python

"""
"""

import os
import pandas as pd
import numpy as np


"""
wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_20/chembl_20_mysql.tar.gz
extract chembl_20_mysql.tar.gz
mysql -p CREATE DATABASE chembl_20;
mysql -p CREATE DATABASE chembl_20;
mysql -u afr -p chembl_20 < chembl_20_mysql/chembl_20.mysqldump.sql

mysql -u afr -p

USE chembl_20;

# get id of compound with a certain smiles
SELECT molregno,canonical_smiles INTO OUTFILE '/tmp/compound_structures.csv' FIELDS TERMINATED BY ',' FROM compound_structures;

# get chembl id of molregno id
SELECT entity_id,chembl_id INTO OUTFILE '/tmp/chembl_id_lookup.csv' FIELDS TERMINATED BY ',' FROM chembl_id_lookup WHERE ENTITY_TYPE = 'compound';

# get chebi ID for compound
SELECT chembl_id,chebi_par_id INTO OUTFILE '/tmp/molecule_dictionary.csv' FIELDS TERMINATED BY ',' FROM molecule_dictionary;
"""

for table in ['compound_structures', 'chembl_id_lookup', 'molecule_dictionary']:
    os.system("cp /tmp/%s.csv ~/Downloads/%s.csv" % (table, table))


df = pd.read_csv("compound_structures.csv", header=None)
df.index = df[1]
df[0].to_csv("compound_structures.indexed.csv")
ref = df.reset_index()
ref.columns = ['smiles', 'id']

df = pd.read_csv("chembl_id_lookup.csv", header=None)
df.index = df[0]
df[1].to_csv("chembl_id_lookup.indexed.csv")
ref1 = df.reset_index()
ref1.columns = ['id', 'chembl']

df = pd.read_csv("molecule_dictionary.csv", header=None)
df = df.replace("\N", np.nan).dropna()
df.index = df[1]
df[0].to_csv("molecule_dictionary.indexed.csv")
ref2 = df.reset_index()
ref2.columns = ['chembl', 'chebi']

ref = pd.merge(pd.merge(ref, ref1, how='outer'), ref2, how='outer')
ref.to_csv("full_mapping.csv", index=False)
ref.dropna().to_csv("full_mapping.existing.csv", index=False)

# Annotate ChEBI entries with function:
os.system("wget ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo")


def get_all_parents_classes(key):
    """
    Get the ChEBI class name/class of every parent of the compound.
    """
    print(key)
    if "is_a" not in ontology[key]:
        return [ontology[key]["name"]]
    else:
        if type(ontology[key]["is_a"]) == list:
            terms = [ontology[key]['name']]
            for term in ontology[key]["is_a"]:
                terms += get_all_parents_classes(term)
            return terms
        elif type(ontology[key]["is_a"]) == str:
            return [ontology[key]['name'], ontology[ontology[key]["is_a"]]["name"]]


def prepare_ontology():
    terms = dict()
    for t in ontology.keys():
        terms[t] = list(set(get_all_parents_classes(t)))


os.system("wget chebi.obo")
ont = parse_obo("chebi.obo")

ontology = dict()
for rec in ont:
    ontology['id'] = rec


# get flatten chebiID:compound_classes mapping
terms = dict()
for t in ontology.keys():
    terms[t] = list(set(get_all_parents_classes(t)))


pd.Series(terms)
