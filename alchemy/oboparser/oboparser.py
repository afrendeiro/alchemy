from collections import defaultdict


__author__ = "Uli Koehler"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__ = "Apache v2.0"


def process_term(goTerm):
    """
    In an object representing a GO term, replace single-element lists with
    their only member.
    Returns the modified object as a dictionary.
    """
    ret = dict(goTerm)  # Input is a defaultdict, might express unexpected behaviour
    for key, value in ret.iteritems():
        if len(value) == 1:
            ret[key] = value[0]
    return ret


def parse_obo(stream):
    """
    Parses a Gene Ontology dump in OBO v1.2 format.
    Yields each
    Keyword arguments:
        stream: The stream to read
    """
    current_term = None
    for line in stream:
        line = line.strip()
        if not line:
            continue  # Skip empty
        if line == "[Term]":
            if current_term:
                yield process_term(current_term)
            current_term = defaultdict(list)
        elif line == "[Typedef]":
            # Skip [Typedef sections]
            current_term = None
        else:  # Not [Term]
            # Only process if we're inside a [Term] environment
            if current_term is None:
                continue
            key, sep, val = line.partition(":")
            current_term[key].append(val.strip())
    # Add last term
    if current_term is not None:
        yield process_term(current_term)
