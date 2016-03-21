import pandas as _pd
from pkg_resources import resource_filename
import os


def load_mappings():
    """
    Load a dataset from alchemy.data

    """
    print resource_filename("alchemy", os.path.join("data", "full_mapping.existing.csv"))
    # return _pd.read_csv(resource_filename("alchemy", os.path.join("data", "full_mapping.existing.csv")))


def load_ontology():
    """
    Load a dataset from alchemy.data

    """
    return _pd.read_csv(resource_filename("alchemy", os.path.join("data", "chebi.obo")))
