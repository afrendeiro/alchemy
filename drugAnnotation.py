#!/usr/bin/env python

"""
Annotate drugs based on SMILE string or drug names from ChemSpider, ChEMBL, ChEBI, KEGG, clinicaltrials.gov,
and create ontology with functional terms from drug annotation.
"""

from argparse import ArgumentParser
from BeautifulSoup import BeautifulStoneSoup
import bioservices
import logging
import numpy as np
import os
import pandas as pd
import re
import sys
import urllib2


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
        description="""Annotate drugs based on SMILE string or drug names from ChemSpider, ChEMBL, ChEBI, KEGG, clinicaltrials.gov,
        and create ontology with functional terms from drug annotation.
        """
    )
    # positional arguments
    parser.add_argument(
        'infile',
        help='CSV file with drug smiles (and optionally names).'
    )
    parser.add_argument(
        'outputFolder',
        help='Folder name to be created and containing output. Will overwrite existing files!'
    )
    # optional arguments
    parser.add_argument(
        '-s', '--smiles-column',
        type=int, default=1, dest='smilesColumn',
        help='Specify the column with the smile strings. Default is column 1 (1-based).'
    )
    parser.add_argument(
        '-n', '--query-name', action='store_true', dest='queryName',
        help='Annotate drugs by name as well as by smile. Off by default.'
    )
    parser.set_defaults(
        queryName=False
    )
    parser.add_argument(
        '--names-column', type=int,  dest='namesColumn',
        help='Specify the column with drug names if available. Default is column 2 (1-based).'
    )
    parser.add_argument(
        '-l', '--logfile', default="log.txt", dest='logFile',
        help='Specify the name of the log file.'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true', dest='verbose',
        help="Verbose behaviour. Will print information to stdout.", default=False
    )

    # parse
    args = parser.parse_args()

    if args.queryName is True and args.namesColumn is None:
        parser.error("--query-name requires you to specify the column number of drug names with --names-column .")
        sys.exit(0)
    elif args.queryName is False and args.namesColumn is not None:
        parser.error("--names-column is only used with the --query-name flag as well.")
        sys.exit(0)

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

    # Make column numbers 0-based
    args.smilesColumn -= 1
    if args.queryName:
        args.namesColumn -= 1

    # Make output directory
    if not os.path.exists(args.outputFolder):
        os.mkdir(args.outputFolder)

    # Infile to pandas data.frame
    if v:
        print("Reading in infile: %s" % args.infile)
    drugInfo = pd.read_csv(args.infile)
    # rename smiles column
    drugInfo.rename(columns={drugInfo.columns[args.smilesColumn]: 'smiles'}, inplace=True)
    # rename names column
    if args.queryName:
        if v:
            print("Going to query by drug name as well as smiles.")
        drugInfo.rename(columns={drugInfo.columns[args.namesColumn]: 'names'}, inplace=True)

    # add drugID column
    drugInfo['drugID'] = ['drug%s' % digit for digit in range(len(drugInfo))]
    drugInfo.to_csv(os.path.join(args.outputFolder, 'drugInfo.csv'), index=False, encoding='utf-8')

    # query ChemSpider
    if v:
        print("Annotating with ChemSpider")
    # csInfo = chemSpiderQuery(drugInfo['smiles'])
    # write to file
    # csInfo.to_csv(os.path.join(args.outputFolder, 'csInfo.csv'), index=False, encoding='utf-8')

    # query ChEMBL
    if v:
        print("Annotating with ChEMBL")
    chemblInfo = chemblQuery(drugInfo['smiles'])
    # write to file
    chemblInfo.to_csv(os.path.join(args.outputFolder, 'chemblInfo.csv'), index=False, encoding='utf-8')

    # query ChEBI
    if v:
        print("Annotating with ChEBI")
    chebiInfo = chebiQuery(drugInfo['smiles'])
    # write to file
    chebiInfo.to_csv(os.path.join(args.outputFolder, 'chebiInfo.csv'), index=False, encoding='utf-8')

    # Annotate ChEMBL info with assays
    if v:
        print("Annotating ChEMBL with assay information")
    chemblAssays = annotateChEMBL(chemblInfo[['drugID', 'chemblID']])
    # write to file
    chemblAssays.to_csv(os.path.join(args.outputFolder, 'chemblAssays.csv'), index=False, encoding='utf-8')

    # query KEGG by name
    if args.queryName:
        if v:
            print("Annotating with KEGG")
        keggInfo = keggQuery(drugInfo['names'])
        # write to file
        keggInfo.to_csv(os.path.join(args.outputFolder, 'keggInfo.csv'), index=False, encoding='utf-8')

    # clinicaltrials.gov annotation
    if args.queryName:
        if v:
            print("Annotating with clinicaltrials.gov")
        ctInfo = queryClinicalTrials(drugInfo['names'])
        ctInfo.to_csv(os.path.join(args.outputFolder, 'ctInfo.csv'), index=False, encoding='utf-8')

    # Concatenate all annotations
    if v:
        print("Concatenating annotations")
    # drugAnnotation = pd.merge(drugInfo, csInfo, how='left', on='drugID')
    drugAnnotation = drugInfo
    drugAnnotation = pd.merge(drugAnnotation, chemblInfo, on='drugID') # , how='left', 
    drugAnnotation = pd.merge(drugAnnotation, chebiInfo, on='drugID') # , how='left', 
    drugAnnotation = pd.merge(drugAnnotation, chemblAssays, on='drugID') # , how='left', 
    if args.queryName:
        drugAnnotation = pd.merge(drugAnnotation, keggInfo, how='left', on='drugID')
        drugAnnotation = pd.merge(drugAnnotation, ctInfo, how='left', on='drugID')

    # save to file
    drugAnnotation.to_csv(os.path.join(args.outputFolder, 'drugAnnotation.csv'), index=False, encoding='utf-8')

    # Create ontology
    if v:
        print("Creating ontology")
    # Remove duplicates based on slimes
    drugAnnotation = drugAnnotation.drop_duplicates('smiles')
    universe = drugAnnotation[pd.notnull(drugAnnotation['chebiOntology'])]
    universe = universe.set_index([range(len(universe))])

    # Export ontology as term-dependent list
    universeSplit = splitOntologyTerms(universe)
    universeSplit.to_csv(os.path.join(args.outputFolder, 'universeOntology.csv'), index=False, encoding='utf-8')

    if v:
        print("Finished.")


def chemSpiderQuery(smiles):
    """
    Queries ChemSpider with SMILES, returns matches.
    Requires list of smiles, outputs pandas DataFrame with hits' information.
    Information is drugID, csID, csName, csSmile.
    """
    from chemspipy import ChemSpider

    CST = os.environ['CHEMSPIDER_SECURITY_TOKEN']
    cs = ChemSpider(security_token=CST)

    csInfo = pd.DataFrame(index=np.arange(len(smiles)), columns=['drugID', 'csID', 'csName', 'csSmile'])

    # request search for each compound based on Smile
    for compound in range(len(smiles)):
        search = cs.async_simple_search(smiles[compound])
        ready = False
        while ready == False:
            if cs.get_async_search_status_and_count(search)['status'] == 'ResultReady':
                ready = True
            elif cs.get_async_search_status_and_count(search)['status'] == 'Failed':
                if v:
                    print("search failed")
                break

        # get tophit
        if cs.get_async_search_result(search) != []:
            tophit = cs.get_async_search_result(search)[0]

            if hasattr(tophit, 'common_name'):
                csInfo.ix[compound] = ['drug%s' % compound, tophit.csid, tophit.common_name, tophit.smiles]
            else:
                csInfo.ix[compound] = ['drug%s' % compound, tophit.csid, np.NaN, tophit.smiles]
        else:
            csInfo.ix[compound] = ['drug%s' % compound] + ['No match'] * 3

    return csInfo


def chemblQuery(smiles):
    """
    Queries ChEMBL with SMILES, returns matches.
    Requires list of smiles, outputs pandas DataFrame with hits' information.
    Information is drugID, chemblID, chemblName, chemblSmile.
    """
    # Initiate ChEMBL request
    chembl = bioservices.ChEMBL()

    # initialize empty DF
    chemblInfo = pd.DataFrame(index=np.arange(len(smiles)), columns=['drugID', 'chemblID', 'chemblName', 'chemblSmile', 'chemblKnown'])

    # check ChEMBL id
    try:
        chemblSmiles = chembl.get_compounds_by_SMILES([smile for smile in smiles])
    except:
        raise

    for i, smile in enumerate(chemblSmiles):
        # check request got match
        if not type(chemblSmiles[i]) == dict:
            if v:
                print('No match for compound %s.' % smile)
            chemblInfo.ix[i] = ['drug%s' % i] + [np.NaN] * 4
        else:
            hits = chemblSmiles[i]['compounds']

            # check there are hits
            if hits != []:
                tophit = hits[0]

                if len(hits) == 1:
                    if v:
                        print("Only one match to a compound found %s. Saved that one." % tophit['chemblId'])
                    chemblInfo.ix[i] = ['drug%s' % i, tophit['chemblId'], tophit['preferredCompoundName'], tophit['smiles'], tophit['knownDrug']]
                else:
                    # If there's more than one match
                    # compare hits (chemical properties)
                    # check if hits have same molecular formula
                    h = []
                    for c in range(len(hits)):
                        # if it's different discard
                        if tophit['molecularFormula'] == hits[c]['molecularFormula']:
                            if v:
                                print("Compound %s has different molecular formula" % hits[c]['chemblId'])
                            h.append(hits[c])
                    hits = h
                    # if only one left, save it
                    if len(hits) == 1:
                        if v:
                            print("Only one match to a compound found %s. Saved that one." % tophit['chemblId'])
                        chemblInfo.ix[i] = ['drug%s' % i, tophit['chemblId'], tophit['preferredCompoundName'], tophit['smiles'], tophit['knownDrug']]
                    else:
                        # check hits have same smile
                        h = []
                        for c in range(len(hits)):
                            # if it's different discard
                            if tophit['smiles'] == hits[c]['smiles']:
                                if v:
                                    print("Compound %s has different smile" % hits[c]['chemblId'])
                                h.append(hits[c])

                        hits = h
                        # if only one left, save it
                        if len(hits) == 1:
                            if v:
                                print("Only one match to a compound found %s. Saved that one." % tophit['chemblId'])
                            chemblInfo.ix[i] = ['drug%s' % i, tophit['chemblId'], tophit['preferredCompoundName'], tophit['smiles'], tophit['knownDrug']]
                        else:
                            # check if there's one which is a known drug
                            h = []
                            for c in range(len(hits)):
                                if hits[c]['knownDrug'] == 'Yes':
                                    if v:
                                        print("Compound %s is known drug" % hits[c]['chemblId'])
                                    h.append(hits[c])

                            hits = h
                            # if not hits remain, add tophit
                            if hits == []:
                                if v:
                                    print("No compound is a known drug. Saving tophit %s." % tophit['chemblId'])
                                chemblInfo.ix[i] = ['drug%s' % i, tophit['chemblId'], tophit['preferredCompoundName'], tophit['smiles'], tophit['knownDrug']]
                            else:
                                if v:
                                    print("%d compounds are known drugs. Saving tophit %s." % (len(hits), tophit['chemblId']))
                                chemblInfo.ix[i] = ['drug%s' % i, tophit['chemblId'], tophit['preferredCompoundName'], tophit['smiles'], tophit['knownDrug']]

    # rename columns to have unique ones
    chemblInfo.rename(columns={'activity': 'chemblActivity', 'targetName': 'targetChemblName', 'organism': 'chemblOrganism'}, inplace=True)
    return chemblInfo


def chebiQuery(smiles):
    """
    Queries ChEBI with SMILES, returns matches.
    Requires list of smiles, outputs pandas DataFrame with hits' information.
    Information is drugID, chebiID, chebiName, chebiSmile.
    """
    import bioservices
    import pandas as pd
    import numpy as np

    # Initiate chebi request
    chebi = bioservices.ChEBI()

    # initialize empty DF
    chebiInfo = pd.DataFrame(index=np.arange(len(smiles)), columns=['drugID', 'chebiID', 'chebiName', 'chebiOntology'])

    for compound, smile in enumerate(smiles):
        # query with smile
        try:
            results = chebi.getStructureSearch(smile, "SMILES", structureSearchCategory = "SIMILARITY", tanimotoCutoff = 0.7)
        except Exception, e:
            if v:
                print("Error. Couldn't query CheBI for compound %s." % smile)
            chebiInfo.ix[compound] = ['drug%s' % compound] + [np.NaN] * 3
            continue

        # check there are results
        if results:
            if v:
                print("Compound %s has results" % smile)
            tophit = results[0][0]
            chebiID = tophit[0]
            chebiName = tophit[1]

            # annotate with ontology
            parentOntology = chebi.getOntologyParents(chebiID)

            # add terms separated by '|'
            terms = ''
            if len(parentOntology) >= 1:
                for parent in parentOntology[0]:
                    if terms == '':
                        terms = str(parent['chebiName'])
                    else:
                        terms = terms + "|" + str(parent[0])

            if v:
                print('Compound drug%s has ontology terms: %s' % (compound, terms))
            # write info
            chebiInfo.ix[compound] = ['drug%s' % compound, chebiID, chebiName, terms]
        else:
            if v:
                print("No results for compound")
            chebiInfo.ix[compound] = ['drug%s' % compound] + [np.NaN] * 3

        #parallel
        #num_cores = multiprocessing.cpu_count()
        #chebiInfo2 = Parallel(n_jobs = num_cores)(delayed(query)(compound, chebiInfo) for compound in xrange(0, len(smiles)))

    return chebiInfo


def keggQuery(names):
    """
    Queries KEGG with drug names, returns matches.
    Requires list of drug strings, outputs pandas DataFrame with hits' information.
    Information is drugID, keggID, keggNames.
    """
    # Initiate ChEMBL request
    kegg = bioservices.KEGG()

    # initialize empty DF
    keggInfo = pd.DataFrame(index=np.arange(len(names)), columns=['drugID', 'keggID', 'keggNames'])

    for compound, name in enumerate(names):
        # if coumpound name is not empty
        if not pd.isnull(names[compound]):
            # get kegg id
            result = kegg.find('drug', names[compound])

            if isinstance(result, unicode) or isinstance(result, str):
                # if there is one result
                if result.strip() != '':
                    keggID = result.strip().split('\t')[0]
                    keggID = re.sub('dr:', '', keggID)
                    keggNames = result.strip().split('\t')[1]
                    keggNames = re.sub(';', '|', keggNames)
                    keggInfo.ix[compound] = ['drug%s' % compound, keggID, keggNames]
                else:
                    keggInfo.ix[compound] = ['drug%s' % compound, np.NaN, np.NaN]
            else:
                keggInfo.ix[compound] = ['drug%s' % compound, np.NaN, np.NaN]
        else:
            keggInfo.ix[compound] = ['drug%s' % compound, np.NaN, np.NaN]

    return keggInfo


def annotateChEMBL(df):
    """
    Annotates matches with ChEMBL ids with function.
    Requires pandas DataFrame with matches from ChEMBL.
    Ouputs same DataFrame fully annotated when possible.
    """
    # Initiate ChEMBL request
    chembl = bioservices.ChEMBL()

    # check ChEMBL id
    chemblActivities = chembl.get_compounds_activities(str(chemblID) for chemblID in df['chemblID'])

    # add annotation columns to DataFrame
    df['activity'] = np.NaN; df['targetChemblId'] = np.NaN; df['targetName'] = np.NaN; df['organism'] = np.NaN

    for i, ID in enumerate(df['chemblID']):
        # check if compound has ChEMBL ID
        if type(ID) is not float or ID != 'No match':
            # if compound has ChEMBL ID, ask if has been used in any assay
            if type(chemblActivities[i]) is dict:
                if chemblActivities[i]['bioactivities'] != []:
                    for assay in chemblActivities[i]['bioactivities']:
                        # check compound is marked as active in assay
                        if assay['activity_comment'] == 'Active':
                            df.loc[i, 'activity'] = 'Active'

                            # check if there is a target ChEMBL ID in this assay
                            if assay['target_chemblid'] != '':
                                df.loc[i, 'targetChemblId'] = assay['target_chemblid']

                            # check if there is a target name in this assay
                            if assay['target_name'] != '' or assay['target_name'] != 'Unspecified' or assay['target_name'] != 'Unchecked':
                                df.loc[i, 'targetName'] = assay['target_name']

                            # check which organism was tested
                            if assay['organism'] != '':
                                df.loc[i, 'organism'] = assay['organism']
                        else:
                            if v:
                                print("Compound %s is not active in assay" % df.ix[i]['drugID'])
                else:
                    if v:
                        print("No assays have been performed with compound %s." % ID)
                    # add placeholder 'No info'
            else:
                if v:
                    print("No assays have been performed with compound %s." % ID)
                # add placeholder 'No info'

        else:
            if v:
                print("Compound %s has no ChEMBL ID." % df['drugID'][i])
            # add placeholder 'No info'

    df.rename(columns={"targetChemblId": "assayTargetChemblId"}, inplace=True)
    df.drop('chemblID', axis=1, inplace=True)

    return df


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


def queryClinicalTrials(names):
    """
    Queries clinicaltrials.gov for trials involving drug.
    Returns data.frame with attributes of trial for completed trials.
    """
    ctInfo = pd.DataFrame(index=np.arange(len(names)), columns=['drugID', 'ctStudyTitle', 'ctStudyUrl', 'ctStudyOutcome'])

    searchURL = 'http://clinicaltrials.gov/search?term='
    studiesURL = 'http://clinicaltrials.gov/show/'

    for compound in range(len(names)):
        if not pd.isnull(names[compound]):
            try:
                url = urllib2.urlopen(searchURL + names[compound] + '&displayxml=true').read()
            except urllib2.HTTPError:
                if v:
                    print('No trials with drug %s.' % names[compound])
                # add NaNs to data.frame
                ctInfo.ix[compound] = ['drug%s' % compound, np.NaN, np.NaN, np.NaN]
                continue

            soup = BeautifulStoneSoup(url)
            if v:
                print("Found %d studies related to drug %s" % (len(soup), names[compound]))

            if soup.findAll('clinical_study') != []:
                # list of strings of studies related to drug
                studies = []; titles = []; urls = []; outcomes = []
                for study in soup.findAll('clinical_study'):
                    # ask if study is not ongoing
                    if study.status.getText() == u'Completed':
                        if v:
                            print('Study %s is completed' % study.findChild('nct_id').getText())
                        # append study ID to list
                        studyID = study.findChild('nct_id').getText()
                        studies.append(studyID)

                        url = urllib2.urlopen(studiesURL + studyID + '?displayxml=true').read()
                        soup = BeautifulStoneSoup(url)

                        if soup.findAll('primary_outcome') != []:
                            titles.append(soup.findAll('brief_title')[0].getText())
                            urls.append(soup.findAll('required_header')[0].findChild('url').getText())
                            outcomes.append(soup.findAll('primary_outcome')[0].findChild('measure').getText())
                        else:
                            if v:
                                print("No conclusive info about study %s" % studyID)

                # add info to data.frame
                ctInfo.ix[compound] = ['drug%s' % compound, '|'.join(titles), '|'.join(urls), '|'.join(outcomes)]

                if v:
                    print("%d studies are completed for drug %s" % (len(studies), names[compound]))
            else:
                # add NaNs to data.frame
                ctInfo.ix[compound] = ['drug%s' % compound, np.NaN, np.NaN, np.NaN]
        else:
            # add NaNs to data.frame
            if v:
                print("No name for compound")
            ctInfo.ix[compound] = ['drug%s' % compound, np.NaN, np.NaN, np.NaN]

    return ctInfo


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
