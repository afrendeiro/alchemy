#!/usr/bin/env python

"""
Test enrichment of terms in a subset of drugs, given annotation of all drugs.
"""

import os, sys, logging, csv, re, bioservices, urllib2
from chemspipy import ChemSpider
import numpy as np
from scipy.stats import fisher_exact
import pandas as pd
import BeautifulSoup as bsoup
from BeautifulSoup import BeautifulStoneSoup

import matplotlib.pyplot as plt
from argparse import ArgumentParser
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
        help = 'CSV file with annotation.')
    # optional arguments
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

    drugAnnotation = pd.io.parsers.read_csv(args.universe)

    # Create ontology
    # Remove duplicates based on slimes
    drugAnnotation = drugAnnotation.drop_duplicates('OUR_SMILES')
    universe = drugAnnotation[pd.notnull(drugAnnotation['chebiOntology'])]
    universe = universe.set_index([range(len(universe))])

    # get file name
    name = re.sub('.csv', '', args.testFile)
    # Read in
    testSet = pd.io.parsers.read_csv(args.testFile)
    # Retrieve only informative columns
    if 'extended' in args.testFile:
        testSet = pd.concat([testSet.ix[:,20], testSet.ix[:,9]], axis=1, keys=['OUR_SMILES', 'score'])
    else:
        testSet = pd.concat([testSet.ix[:,2], testSet.ix[:,12], testSet.ix[:,13]], axis=1, keys=['OUR_SMILES', 'change', 'p-value'])
    # Remove duplicates based on slimes
    testSet = testSet.drop_duplicates('OUR_SMILES')
    # Remove rows without smile
    testSet = testSet[pd.notnull(testSet['OUR_SMILES'])]
    # annotate testSet files with drugAnnotation
    testSet = pd.merge(testSet, drugAnnotation, how = 'left', on = 'OUR_SMILES')
    # exclude drugs not annotated
    testSet = testSet[pd.notnull(testSet['chebiOntology'])]

    # filter by p-value (optional)
    #testSetSetPos = testSet[testSet['p-value'] < 0.01]
    if 'extended' in args.testFile:
        testSet = testSet[['cemmID', 'chebiID', 'chebiOntology', 'DRUG_NAME']]
        print("Drugs enriched in file %s:" % args.testFile)
        testSetEnrich = testEnrichment(universe, testSet)

        # Export significant data
        testSetEnrichExport = pd.DataFrame(testSetEnrich).T
        testSetEnrichExport.columns = ['odds_ratio', 'p-value', 'counts', 'drugs with term']
        testSetEnrichExport = testSetEnrichExport.sort('p-value')
        testSetEnrichExport.to_csv(name + '_tested.csv', index = True, encoding='utf-8')
    else:
        # separate positive and negative changes
        testSetSetPos = testSet[testSet['change'] > 0]
        testSetSetPos = testSetSetPos[['cemmID', 'chebiID', 'chebiOntology', 'DRUG_NAME']]
        testSetSetNeg = testSet[testSet['change'] < 0]
        testSetSetNeg = testSetSetNeg[['cemmID', 'chebiID', 'chebiOntology', 'DRUG_NAME']]

        print("Positive drugs enriched in file %s:" % args.testFile)
        testSetPosEnrich = testEnrichment(universe, testSetSetPos)
        print("Negative drugs enriched in file %s:" % args.testFile)
        testSetNegEnrich = testEnrichment(universe, testSetSetNeg)

        # Export significant data
        testSetEnrichPosExport = pd.DataFrame(testSetPosEnrich).T
        testSetEnrichPosExport['direction'] = 'positive'
        testSetEnrichNegExport = pd.DataFrame(testSetNegEnrich).T
        testSetEnrichNegExport['direction'] = 'negative'

        testSetEnrichExport = pd.concat([testSetEnrichPosExport, testSetEnrichNegExport])
        testSetEnrichExport.columns = ['odds_ratio', 'p-value', 'counts', 'drugs with term', 'direction']
        testSetEnrichExport = testSetEnrichExport.sort('p-value')
        testSetEnrichExport.to_csv(name + '_tested.csv', index = True, encoding='utf-8')


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
                SET.append((term, df['chebiID'][i], df['cemmID'][i]))

    SET = pd.DataFrame(SET)
    SET.columns = ['term', 'chebiID', 'cemmID']
    return SET
    
    
def countDrugsWithTerm(drugSet, term):
    """
    DEPRECATED!.
    """
    counts = [0,0]
    for drug in drugSet['term']:
        if drug == term:
            counts[0] += 1
        else:
            counts[1] += 1
    return counts


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
                drugs.append((drugSet['cemmID'][drug], drugSet['DRUG_NAME'][drug]))
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