#!/usr/bin/env python

"""
This is the main code for the alchemy package.
"""

import sys
import pandas as pd
from argparse import ArgumentParser
from oboparser import parse_obo


def uniquify(function):
    """
    Decorator to make the output of functions returning lists, unique.
    """
    def wrapper(obj, *args):
        return list(set(function(obj, *args)))
    return wrapper


class IdMapper(object):
    """
    Map smiles to compound IDs.
    """
    def __init__(self):
        self.read_dbs()

    def read_dbs(self):
        """
        Read representation of ChEMBL database and prepare ID mappings.
        """
        import os
        from pkg_resources import resource_stream
        self.mapping = pd.read_csv(resource_stream('alchemy', os.path.join('data', "full_mapping.existing.csv")))

    def annotate(self, keys, id_type="smiles"):
        types = ["smiles", "id", "chembl", "chebi", "name"]
        if id_type in types:
            return pd.merge(keys, self.mapping, how="left")
        else:
            raise TypeError("Provided type of input is not supported. Please choose on of: %s" % ",".join(types))
        # self.chebi_ids = Parallel(n_jobs=self.parameters.processors)(delayed(self.match_smile)(smile) for smile in smiles)


class Ontology(object):
    """
    Ontology handler for alchemy.
    """
    def __init__(self):
        self.parse_obo()

    def __repr__(self):
        return "Ontology handler."

    def parse_obo(self):
        """
        Read representation of ChEMBL database and prepare ID mappings.
        """
        from pkg_resources import resource_stream

        obo = parse_obo(resource_stream("alchemy", "data/chebi.obo"))

        self.ontology = dict()
        for entry in obo:
            self.ontology[entry["id"]] = entry

    @uniquify
    def get_all_parents_classes(self, key):
        """
        Get the ChEBI class name/class of every parent of the compound.
        """
        if key not in self.ontology:
            print("Warning. You're querying non-existant ChEBI ids! %s" % key)
            return []
        if "is_a" not in self.ontology[key]:
            return [self.ontology[key]["name"]]
        else:
            if type(self.ontology[key]["is_a"]) == list:
                terms = [self.ontology[key]["name"]]
                for term in self.ontology[key]["is_a"]:
                    terms += self.get_all_parents_classes(term)
                return terms
            elif type(self.ontology[key]["is_a"]) == str:
                return [self.ontology[key]["name"]] + self.get_all_parents_classes(self.ontology[key]["is_a"])

    @uniquify
    def get_all_parent_roles(self, key):
        """
        Get the ChEBI function of a compound and of all of its parent functions.
        """
        import re
        if key not in self.ontology:
            print("Warning. You're querying non-existant ChEBI ids! %s" % key)
            return []
        if "relationship" not in self.ontology[key]:
            return []
        else:
            if type(self.ontology[key]["relationship"]) == list:
                terms = []
                for term in self.ontology[key]["relationship"]:
                    if "has_role " in term:
                        terms += self.get_all_parent_roles(re.sub("has_role ", "", term))
                return terms
            elif type(self.ontology[key]["relationship"]) == str:
                term = self.ontology[key]["relationship"]
                if "has_role " in term:
                    return self.get_all_parent_roles(re.sub("has_role ", "", term))
                else:
                    return []

    def get_all_classes(self):
        """
        Get the ChEBI class name/class of every parent of the compound.
        """
        all_terms = list()
        for term in self.ontology.keys():
            all_terms += self.get_all_parents_classes(term)
        return all_terms

    def get_all_roles(self):
        """
        Get any known roles of every parent of the compound.
        """
        all_terms = list()
        for term in self.ontology.keys():
            all_terms += self.get_all_parent_roles(term)
        return all_terms


def add_args(parser):
    """
    """
    parser.add_argument(dest="input_file", type=str,
                        help="Input file with query.")
    parser.add_argument("-b", "--background", required=False,
                        dest="background", type=str,
                        help="Background set to compare enrichment to. Must be in same format as input file.")
    parser.add_argument("-t", "--input-type",
                        dest="input_type", type=str, default="smiles",
                        choices=["smiles", "id", "chembl", "chebi", "name"],
                        help="Type of input query provided.")
    parser.add_argument("-d", "--delimiter",
                        dest="delimiter", type=str, default=",",
                        help="Delimiter character in input file.")
    parser.add_argument("--header",
                        dest="header", action="store_true",
                        help="Does the input file contain a header?")
    parser.add_argument("-a", "--annotation-output", dest="annotation_output", default="alchemy_compound_annotation.csv",
                        help="Directory to write results to. Will create if does not exist, will overwrite if existing.")
    parser.add_argument("-ab", "--background-annotation-output", dest="background_annotation_output", default="alchemy_background_annotation.csv",
                        help="Directory to write results to. Will create if does not exist, will overwrite if existing.")
    parser.add_argument("-e", "--enrichment-ouput", dest="enrichment_output", default="alchemy_compound_class_enrichment.csv",
                        help="Directory to write results to. Will create if does not exist, will overwrite if existing.")
    parser.add_argument("--no-enrich", dest="enrich", action="store_false", default=True,
                        help="If enrichment analysis is to be performed.")
    parser.add_argument("-p", "--processors", dest="processors", type=int,
                        help="Number of parallel processors to use.")
    parser.add_argument("--tmp-dir", dest="tmp_dir", type=str, default='/tmp/joblib',
                        help="Temporary directory to write cache to.")

    return parser


def parse_input(args, input_attr="input_file"):
    df = pd.read_table(getattr(args, input_attr), sep=args.delimiter, header=0 if args.header else None)
    df.columns = [args.input_type]
    return df


def number_to_chebi(chebi_number):
    try:
        return "CHEBI:%s" % str(int(chebi_number))
    except ValueError:
        return pd.np.nan


args = "-t name -a annotation.csv -e enrichment.csv ~/drugs.txt".split(" ")


def main(args=None):
    # Parse arguments
    parser = ArgumentParser(
        prog=__name__,
        description="%s. Annotate your compounds and explore enriched classes/functional roles." % __name__
    )
    parser = add_args(parser)

    # Parse
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    # Parse user input
    keys = parse_input(args)

    # Parse ChEMBL dbs
    ids = IdMapper()

    # Get ChEBI ids
    annotated_query = ids.annotate(keys, id_type=args.input_type)
    # bring along keys without match to database
    annotated_query['chebi'] = map(number_to_chebi, annotated_query['chebi'])

    if annotated_query.dropna().shape[0] < 1:
        raise KeyError("No provided compound could be matched to database.")

    # Parse ChEBI ontology
    ontology = Ontology()

    # Get the chemical classes of the parents of each smiles
    annotated_query["chebi_classes_parents"] = annotated_query.apply(
        lambda x: ",".join(ontology.get_all_parents_classes(x['chebi'])) if x['chebi'] is not pd.np.nan else pd.np.nan,
        axis=1
    )
    # Get any functions of the parents of each smiles
    annotated_query["chebi_roles_parents"] = annotated_query.apply(
        lambda x: ",".join(ontology.get_all_parent_roles(x['chebi'])) if x['chebi'] is not pd.np.nan else pd.np.nan,
        axis=1
    )

    # Annotate with name
    annotated_query["name"] = annotated_query.apply(
        lambda x: ontology.ontology[x['chebi']]['name'] if x['chebi'] is not pd.np.nan else pd.np.nan,
        axis=1
    )

    # Annotate with ids of other references
    annotated_query["references"] = annotated_query.apply(
        lambda x: ontology.ontology[x['chebi']]['xref'] if x['chebi'] is not pd.np.nan else pd.np.nan,
        axis=1
    )

    # Output annotation
    annotated_query.to_csv(args.annotation_output, index=False)

    # End here if enrichment test is not required
    if not args.enrich:
        return

    # Make long pandas dataframe with:
    # columns:
    # ["smiles", "chebi_id", "classes_roles"]
    c = annotated_query['chebi_classes_parents'].str.split(',').apply(pd.Series, 1).stack()
    c.index = c.index.droplevel(-1)
    c.name = "classes_roles"
    classes = annotated_query.join(c).dropna()

    r = annotated_query['chebi_roles_parents'].str.split(',').apply(pd.Series, 1).stack()
    r.index = r.index.droplevel(-1)
    r.name = "classes_roles"
    roles = annotated_query.join(r).dropna()

    # "classes_roles" now represents both things in equal terms
    # this should be further understood and perhaps filtered for some common terms
    # and tested independently
    query_classes = classes.append(roles)

    # count number of compounds per class
    query_count = query_classes.groupby("classes_roles").apply(len)

    # Select background/universe:
    if args.background is None:
        # if not provided, use whole database
        # warn about p-value
        from collections import Counter
        univ_count = pd.Series(Counter(ontology.get_all_classes())).append(pd.Series(Counter(ontology.get_all_roles())))
    else:
        # if provided, use that - report amount of overlap
        # Parse user input
        keys = parse_input(args, "background")

        # TODO:
        # report amount of overlap between query and background
        # warn if query not in background

        # Get ChEBI ids
        annotated_background = ids.annotate(keys)
        # bring along keys without match to database
        annotated_background['chebi'] = map(number_to_chebi, annotated_background['chebi'])

        # Get the chemical classes of the parents of each smiles
        annotated_background["chebi_classes_parents"] = annotated_background.apply(
            lambda x: ",".join(ontology.get_all_parents_classes(x['chebi'])) if x['chebi'] is not pd.np.nan else pd.np.nan,
            axis=1
        )
        # Get any roles of the parents of each smiles
        annotated_background["chebi_roles_parents"] = annotated_background.apply(
            lambda x: ",".join(ontology.get_all_parent_roles(x['chebi'])) if x['chebi'] is not pd.np.nan else pd.np.nan,
            axis=1
        )

        # write annotated output
        annotated_background.to_csv(args.background_annotation_output, index=False)

        # Make long pandas dataframe with:
        # columns:
        # ["smiles", "chebi_id", "classes"]
        uc = annotated_background['chebi_classes_parents'].str.split(',').apply(pd.Series, 1).stack()
        uc.index = uc.index.droplevel(-1)
        uc.name = "classes_roles"
        univ_classes = annotated_background.join(uc).dropna()

        ur = annotated_background['chebi_roles_parents'].str.split(',').apply(pd.Series, 1).stack()
        ur.index = ur.index.droplevel(-1)
        ur.name = "classes_roles"
        univ_roles = annotated_background.join(ur).dropna()

        # "classes_roles" now represents both things in equal terms
        # this should be further understood and perhaps filtered for some common terms
        # and tested independently
        univ = univ_classes.append(univ_roles)
        # count number of compounds per class
        univ_count = univ.groupby("classes_roles").apply(len)

    # build contingency tables
    counts = pd.concat([query_count, univ_count], axis=1)
    counts.columns = ["query", "universe"]

    # drop nan to not test classes without compounds in it
    counts = counts.dropna()

    # total of each group
    f_sum = counts["query"].sum()
    u_sum = counts["universe"].sum()

    contigency_table = counts.apply(lambda x: pd.Series([
        x["query"],
        x["universe"],
        f_sum - x["query"],
        u_sum - x["universe"]]), axis=1)
    contigency_table.columns = ["support", "b", "c", "d"]

    # test overrepresentation
    from scipy.stats import fisher_exact
    stats_table = contigency_table.apply(lambda x: pd.Series(fisher_exact([[x["support"], x["b"]], [x["c"], x["d"]]])), axis=1)
    stats_table.columns = ["oddsratio", "pvalue"]

    # correct pvalues
    from statsmodels.sandbox.stats.multicomp import multipletests
    stats_table["corrected_pvalue"] = multipletests(stats_table["pvalue"])[1]

    # output csv
    output = contigency_table.join(stats_table)
    try:
        output.sort_values(by=['corrected_pvalue'], inplace=True)
    except:
        output.sort(['corrected_pvalue'], inplace=True)

    output.to_csv(args.enrichment_output, index=True)

    if __name__ != "__main__":
        return (annotated_query, output)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
