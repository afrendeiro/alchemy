#!/usr/bin/env python

"""
Test enrichment of terms in a subset of drugs, given annotation of all drugs.
"""

from argparse import ArgumentParser
import os, sys, logging, csv, re
import numpy as np
from scipy.stats import fisher_exact
import pandas as pd
from collections import Counter

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2014, Andre F. Rendeiro"
__credits__ = []
__license__ = "GPL3"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"

def main():
    # argparser    
    parser = ArgumentParser(
        description = 'Test enrichment of terms in a subset of drugs, given annotation of all drugs.',
        usage       = 'python drugEnrichment.py [OPTIONS] testFile.csv universeFile.csv'
        )
    # positional arguments
    parser.add_argument('testFile',
        help = 'CSV file with drugs to test.')
    parser.add_argument('universe',
        help = 'CSV file with annotation of all drugs (universe).')
    # optional arguments
    parser.add_argument('-t', '--id-column-test', type = int, default = 1, dest = 'idColumnTest',
        help = 'Specify the column with an ID string for the test file. Default is column 1 (1-based).')
    parser.add_argument('-u', '--id-column-univ', type = int, default = 1, dest = 'idColumnUniv',
        help = 'Specify the column with an ID string for the universe file. Default is column 1 (1-based).')
    parser.add_argument('-o', '--ontology-column', type = int, default = 3, dest = 'ontologyColumn',
        help = 'Specify the column with the ontology strings in the universe file. Default is column 3 (1-based).')
    parser.add_argument('-n', '--ontology-name', action='store_true', dest = 'names',
        help = 'Ontology file has name column. Will be used to when displaying output. Off by default.')
    parser.set_defaults(queryName = False)
    parser.add_argument('-s', '--score-column', type = int, default = 0, dest = 'scoreColumn',
        help = 'OPTIONAL: Specify a column with a score to test enrichment dependent on signal of the score (will split entries and perform test seperately for entries with positive and negative score) - 1-based.')
    parser.add_argument('-l', '--logfile', default = "log.txt", dest = 'logFile',
        help = 'Specify the name of the log file.')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose',
        help = "Verbose behaviour. Will print information to stdout.", default = False)
    # parse
    args = parser.parse_args()

    # logging
    global logger
    logger = logging.getLogger(sys.argv[0])
    logger.setLevel(logging.INFO)

    # create a file handler
    handler = logging.FileHandler(args.logFile)
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)

    global v
    v = args.verbose

    global names
    names = args.names

    # Make column numbers 0-based
    args.idColumnTest -= 1
    args.idColumnUniv -= 1
    args.ontologyColumn -= 1

    universe = pd.io.parsers.read_csv(args.universe)

    ### Create ontology
    # rename drugID column
    universe.rename(columns = {universe.columns[args.idColumnUniv]: 'drugID'}, inplace = True)
    # Remove duplicates based on drugID
    universe = universe.drop_duplicates('drugID')
    # rename ontology column
    universe.rename(columns = {universe.columns[args.ontologyColumn]: 'chebiOntology'}, inplace = True)
    universe = universe[pd.notnull(universe['chebiOntology'])]
    universe = universe.set_index([range(len(universe))])

    ### Test file
    testName = re.sub('.csv', '', args.testFile)
    testSet = pd.io.parsers.read_csv(args.testFile)
    # rename drugID column
    testSet.rename(columns = {testSet.columns[args.idColumnTest]: 'drugID'}, inplace = True)
    # Remove duplicates based on slimes
    testSet = testSet.drop_duplicates('drugID')
    # Remove rows without drugID
    testSet = testSet[pd.notnull(testSet['drugID'])]
    # annotate testSet files with universe
    if names:
        testSet = testSet[['drugID', 'names']]
    else:
        testSet = pd.DataFrame(testSet['drugID'])

    testSet = pd.merge(testSet, universe, how = 'left', on = 'drugID')
    
    # exclude drugs not annotated
    testSet = testSet[pd.notnull(testSet['chebiOntology'])]

    if args.scoreColumn == 0:
        if v:
            print("Drugs enriched in file %s:" % args.testFile)
        testSetEnrich = testEnrichment(universe, testSet)

        # Export tested data
        testSetEnrichExport = pd.DataFrame(testSetEnrich).T
        testSetEnrichExport.columns = ['odds_ratio', 'p-value', 'counts', 'drugs with term']
        testSetEnrichExport = testSetEnrichExport.sort('p-value')
        testSetEnrichExport.to_csv(testName + '_tested.csv', index = True, encoding='utf-8')
    else:
        # rename score column
        testSet.rename(columns = {testSet.columns[args.scoreColumn - 1]: 'score'}, inplace = True)

        # separate positive and negative scores
        testSetSetPos = testSet[testSet['score'] > 0]
        testSetSetNeg = testSet[testSet['score'] < 0]
        if v:
            print("Positive drugs enriched in file %s:" % args.testFile)
        testSetPosEnrich = testEnrichment(universe, testSetSetPos)
        if v:
            print("Negative drugs enriched in file %s:" % args.testFile)
        testSetNegEnrich = testEnrichment(universe, testSetSetNeg)

        # Export tested data
        testSetEnrichPosExport = pd.DataFrame(testSetPosEnrich).T
        testSetEnrichPosExport['direction'] = 'positive'
        testSetEnrichNegExport = pd.DataFrame(testSetNegEnrich).T
        testSetEnrichNegExport['direction'] = 'negative'

        testSetEnrichExport = pd.concat([testSetEnrichPosExport, testSetEnrichNegExport])
        testSetEnrichExport.columns = ['odds_ratio', 'p-value', 'counts', 'drugs with term', 'direction']
        testSetEnrichExport = testSetEnrichExport.sort('p-value')
        testSetEnrichExport.to_csv(testName + '_tested.csv', index = True, encoding='utf-8')


def splitOntologyTerms(df):
    """
    Split pandas dataframe chebiontology terms and returns dataframe with ids and terms.
    Usefull to export.
    """
    # Export ontology as term-dependent list
    df = df.set_index([range(len(df))])
    SET = []
    for i in range(len(df.chebiOntology)):
        string = df.chebiOntology[i]
        if not pd.isnull(string):
            terms = string.split('|')
            for term in terms:
                SET.append((df['drugID'][i], df['chebiID'][i], term))

    SET = pd.DataFrame(SET)
    SET.columns = ['drugID', 'chebiID', 'term']
    return SET
    

def testEnrichment(universe, drugSet):
    """
    Counts number of drugs annotated with term in set and tests significance with Fisher's exact test.
    """
    drugSet = drugSet.set_index([range(len(drugSet))])
    # random sample from universe
    #drugSet = universe.ix[random.sample(universe.index, 100)]

    # get unique terms in the universe
    termSet = splitOntologyTerms(drugSet).term.unique()

    # count occurrences in subset and universe
    termCounts = {}
    for term in termSet:
        counts = [[0,0],[0,0]]
        drugs = []

        for drug in range(len(drugSet)):
            if term in drugSet.chebiOntology[drug].split('|'):
                counts[0][0] += 1
                if names:
                    drugs.append((drugSet['drugID'][drug], drugSet['names'][drug]))
                else:
                    drugs.append(drugSet['drugID'][drug])
            else:
                counts[0][1] += 1

        for drug in range(len(universe)):
            if term in universe.chebiOntology[drug].split('|'):
                counts[1][0] += 1
            else:
                counts[1][1] += 1
        
        termCounts[term] = (counts, drugs)

    # test significance
    testedTerms = {}
    for term, (counts, drugs) in termCounts.iteritems():
        oddsratio, pvalue = fisher_exact(counts)
        if v:
            print("Term '%s' is significantly enriched. P-value: %f" % (term, pvalue), counts, "Drugs in set: %s" % drugs)
        testedTerms[term] = [oddsratio, pvalue, counts, drugs]

    return dict(sorted(testedTerms.items(), key=lambda e: e[1][1]))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)
