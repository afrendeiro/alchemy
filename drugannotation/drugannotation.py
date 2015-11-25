#!/usr/bin/env python

"""
This is the main script of the cll-patients project.
"""

import sys
import pandas as pd
from argparse import ArgumentParser
from joblib import Memory, Parallel, delayed
import obo


class InputParser(object):
    """
    Parse file with SMILES into entries.
    """
    def __init__(self, args):
        self.smiles = pd.read_table(args.input_file, sep=self.args.delimiter)

    def __call__(self):
        return self.smiles


class IdMapper(object):
    """
    """
    def __init__(self, args):
        self.args = args
        self.read_dbs()

    def read_dbs(self):
        """
        Read representation of ChEMBL database and prepare ID mappings.
        """
        self.id_chembl_mapping = 1
        self.id_chembl_mapping = 1
        """
        # get id of compound with a certain smiles
        SELECT MOLREGNO FROM COMPOUND_STRUCTURES WHERE CANONICAL_SMILES = 'smiles'

        # get chembl id of molregno id
        SELECT CHEMBL_ID FROM CHEMBL_ID_LOOKUP WHERE ENTITY_ID = 'molregno' AND ENTITY_TYPE = 'compound'

        # get chebi ID for compound
        SELECT CHEBI_PAR_ID FROM MOLECULE_DICTIONARY WHERE CHEMBL_ID = 'chemblid'
        """

    def match_smiles(self, smiles):
        self.chebi_ids = Parallel(n_jobs=self.parameters.processors)(delayed(sum)(i ** 2) for i in range(10))

        return self.chebi_ids


class Ontology(object):
    """
    """
    def __init__(self):
        pass

    def parse_obo(self, ontology_file):
        """
        Read representation of ChEMBL database and prepare ID mappings.
        """
        parser = obo.Parser(open(ontology_file, "r"))
        self.ontology = dict()
        for stanza in parser:
            self.ontology[stanza.tags["id"][0]] = stanza.tags

    def get_all_parents_classes(self, key, parents=list()):
        """
        Get the ChEBI class name/class of every parent of the compound.
        """
        if not hasattr(self.ontology[key][1], "is_a"):
            return [self.ontology[k][1]["name"] for k in parents + [key]]  # get parents and the term itself
        else:
            self.get_all_parents(key, parents + self.ontology[key][1]["is_a"])


def main():
    # Parse arguments
    parser = ArgumentParser(
        prog="pipelines",
        description="pipelines. Project management and sample loop."
    )
    parser = add_args(parser)

    # Parse
    args = parser.parse_args()

    # tmp joblib dir
    mem = Memory(cachedir=args.tmp_dir)

    # Parse user input
    smiles = InputParser(args)()

    # Parse ChEMBL dbs
    ids = IdMapper()

    # Get ChEBI ids
    chebi_ids = ids.match_smiles(smiles)

    # Parse ChEBI ontology
    ontology = Ontology()
    ontology.parse_obo("")

    # Get the chemical classes of the parents of each smiles
    chebi_parents = {key: ontology.get_all_parents_classes(key) for key in chebi_ids}

    # Make long pandas dataframe with:
    # columns:
    # ["smiles", "chebi_id", "classes"]

    # groupby "classes"

    # Select background/universe:
    # if provided, use that - report amount of overlap

    # if not use whole database?

    # build confusion matrices
    # test overrepresentation

    # output csv


def add_args(parser):
    """
    """
    parser.add_argument(dest="input_file", type=str,
                        help="Input file with SMILES.")
    parser.add_argument("-d", "--delimiter",
                        dest="delimiter", type=str, default=",",
                        help="Delimiter character in input file.")
    parser.add_argument("-o", "--output-dir", dest="output_dir",
                        help="Directory to write results to. Will create if does not exist, will overwrite if existing.")
    parser.add_argument("-p", "--processors", dest="processors", type=int,
                        help="Number of parallel processors to use.")
    parser.add_argument("--tmp-dir", dest="tmp_dir", type=int, default='/tmp/joblib',
                        help="Temporary directory to write cache to.")

    return parser


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
