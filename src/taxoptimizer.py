#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# corinne.maufrais@pasteur.fr
#

# version 2.0


import os
import sys
import argparse
from bsddb import db as bdb

import Golden

try:
    GOLDENDATA = os.environ['GOLDENDATA']
except:
    GOLDENDATA = "/local/gensoft2/exe/golden/1.1a/share/golden/db/"
    # GOLDENDATA = "/mount/banques/prod/index/golden/"
    os.environ['GOLDENDATA'] = GOLDENDATA


# ###################


class ParserError:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[ParserError] " + self.err


class ParserWarning:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[ParserWarning] " + self.err


class InputLine:
    # vlegrand@pasteur.fr development
    def __init__(self, orig_line, acc, skip_db):
        self.orig_line = orig_line
        self.acc = acc
        self.skip_db = skip_db


def parseUniprot(flatFile, DE):
    """
    parse uniprot or embl like flat file
    """
    description = ''
    taxId = ''
    taxoLight = ''
    orgName = ''
    lineFld = flatFile.split('\n')
    vuOS = False
    vuOC = False
    vuOCXX = False
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:2].strip()
        value = line[5:]
        value = value.replace('\n', '')
        if DE and tag == 'DE':
            description += value
        elif tag == 'OS' and not vuOS:
            fldOS = value.split('(')
            orgName += fldOS[0].strip()
            vuOS = True
        elif tag == 'OC' and not vuOCXX:
            taxoLight += value
            vuOC = True
        elif tag == 'XX' and vuOC:
            vuOCXX = True
        elif tag == 'OX':
            taxId += value
        elif tag in ['RN', 'DR', 'CC', 'FH', 'SQ']:
            return orgName, taxId, taxoLight, description
    return orgName, taxId, taxoLight, description


def parseGenbank(flatFile, DE):
    """
    parse genbank like flat file
    """
    description = ''
    taxId = ''
    taxoLight = ''
    orgName = ''
    taxoL = False
    lineFld = flatFile.split('\n')
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:12].strip()
        value = line[12:]
        value = value.replace('\n', '')
        if DE and tag == 'DEFINITION':
            description += value
        elif tag == 'ORGANISM':  # one line ##### pb qqfois plusieurs ligne
            orgName += value.strip()
            taxoL = True
        elif taxoL and not tag:  # multiple line without tag
            taxoLight += value
        elif tag == 'TaxID':  # gi number
            taxId += value
        elif tag in ['REFERENCE', 'COMMENT', 'FEATURES', 'ORIGIN']:
            return orgName, taxId, taxoLight, description
        elif not tag:  # fin taxonomy
            taxoL = False
    return orgName, taxId, taxoLight, description


def parse(flatFile, DE):
    """
    parse db flat file (uniprot, embl, genbank) into a DBRecord object
    """
    if flatFile[:2] == 'ID':
        return parseUniprot(flatFile, DE)
    elif flatFile[:5] == 'LOCUS':
        return parseGenbank(flatFile, DE)
    return '', '', '', ''


class TaxOptimizerError:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[taxoptimizer] " + self.err

##############################################################################
#
#            Golden: Keep this for compatibility and performance testing.
#
##############################################################################


def doGolden(db, ac, DE):
    # ########################## db ref
    if db in['sp', 'sw', 'swissprot', 'tr', 'trembl']:
        db = 'uniprot'
    elif db in ['emb', 'dbj']:
        db = 'embl'
    elif db in ['gb']:
        db = 'genbank'
    elif db in ['gp']:
        db = 'genpept'
    elif db in ['ref']:
        db = 'refseq'
    elif db in ['rdpii']:
        db = 'rdpii'
    elif db[0:8] == 'embl_wgs':
        db = 'embl_wgs'
    elif db[0:8] == 'genbank_wgs':
        db = 'genbank_wgs'
    elif db in ['pir', 'pdb', 'tpg', 'tpe', 'tpd', 'prf']:
        return '', '', '', ''
    try:
        flatFile = Golden.access(db, ac)
    except IOError, err:
        print >>sys.stderr, err, db, ac
        sys.exit()
    if flatFile:
        orgName, taxId, taxoLight, description = parse(flatFile, DE)  # orgName, taxId, taxoLight, description
        flatFile = ''  # buffer free
        return orgName, taxId, taxoLight, description
    else:
        return '', '', '', ''


# Builds query input string
def buildQueryStr(txt_line, db, acc, l_cards, cnt_cards, l_lines):
    # vlegrand@pasteur.fr development
    skip_db = False
    if db in['sp', 'sw', 'swissprot', 'tr', 'trembl']:
        db = 'uniprot'
    elif db in ['emb', 'dbj']:
        db = 'embl'
    elif db in ['gb']:
        db = 'genbank'
    elif db in ['gp']:
        db = 'genpept'
    elif db in ['ref']:
        db = 'refseq'
    elif db in ['rdpii']:
        db = 'rdpii'
    elif db[0:8] == 'embl_wgs':
        db = 'embl_wgs'
    elif db[0:8] == 'genbank_wgs':
        db = 'genbank_wgs'
    elif db in ['pir', 'pdb', 'tpg', 'tpe', 'tpd', 'prf']:
        # return '', '', '', ''
        skip_db = True
    if not skip_db:
        l_cards += db
        l_cards += ":"
        l_cards += acc
        l_cards += "\n"
    cnt_cards += 1
    li = InputLine(txt_line, acc, skip_db)
    l_lines.append(li)
    return l_cards, cnt_cards, l_lines

##############################################################################
#
#            Golden Multi. Using new version of golden.
#            #vlegrand@pasteur.fr
#
#  l_input : Depending on input_flag, a list of bank:AC bank:locus separated by '\n'
#  so, bank:AC\nbank:locus...
##############################################################################


def doGoldenMulti(allTaxo, l_cards, DE, allTaxId, ncbi_osVSocBDB):
    idx_res = 0
    lst_input = l_cards.split("\n")
    try:
        flatFile = Golden.access_new(l_cards)
        while (flatFile is not None):
            db_acc = lst_input[idx_res]
            l_db_acc = db_acc.split(":")
            acc = l_db_acc[1]
            if acc not in allTaxo:
                allTaxo[acc] = {'db': db}
                allTaxo[acc]['orgName'], allTaxo[acc]['taxId'], allTaxo[acc]['taxoLight'], allTaxo[acc]['DE'] = parse(flatFile, DE)  # orgName, taxId, taxoLight, description
                allTaxo, allTaxId = extractTaxoFrom_osVSocBDB_multi(acc, allTaxo, allTaxId, ncbi_osVSocBDB)

            flatFile = Golden.access_new(l_cards)
            idx_res += 1
        return allTaxo
    except IOError, err:
        print >>sys.stderr, err, l_cards
        sys.exit()


# display results from allTaxo dictionnary.
def printResults(l_lines, allTaxo, outfh, notaxfhout, splitFile):
    for li in l_lines:
        taxonomy = ''
        if not li.skip_db:
            if 'taxoFull' in allTaxo[li.acc]:
                taxonomy = allTaxo[li.acc]['taxoFull']
            else:
                taxonomy = allTaxo[li.acc]['taxoLight']

        if taxonomy:
            print >>outfh, li.orig_line, "\t%s\t%s\t%s" % (allTaxo[li.acc]['orgName'], taxonomy, allTaxo[li.acc]['DE'])
        else:
            if notaxfhout:
                print >>notaxfhout, li.orig_line
            if not splitFile:
                print >>outfh, li.orig_line

##############################################################################
#
#            Taxonomy
#
##############################################################################


def extractTaxoFrom_osVSocBDB(acc, allTaxo, allTaxId, BDB):
    taxonomy = allTaxo[acc]['taxoLight']
    orgName = allTaxo[acc]['orgName']
    if orgName and orgName not in allTaxId:
        taxoFull = BDB.get(str(orgName))
        if taxoFull:
            allTaxo[acc]['taxoFull'] = taxoFull
            allTaxId[orgName] = taxoFull
            allTaxo[acc]['taxoLight'] = ''
            taxonomy = taxoFull
    elif orgName:
        allTaxo[acc]['taxoFull'] = allTaxId[orgName]
        allTaxo[acc]['taxoLight'] = ''
        taxonomy = allTaxId[orgName]
    return taxonomy, allTaxo, allTaxId


def extractTaxoFrom_osVSocBDB_multi(acc, allTaxo, allTaxId, BDB):
    orgName = allTaxo[acc]['orgName']
    if orgName and orgName not in allTaxId:
        taxoFull = BDB.get(str(orgName))
        if taxoFull:
            allTaxo[acc]['taxoFull'] = taxoFull
            allTaxId[orgName] = taxoFull
            allTaxo[acc]['taxoLight'] = ''
    elif orgName:
        allTaxo[acc]['taxoFull'] = allTaxId[orgName]
        allTaxo[acc]['taxoLight'] = ''
    return allTaxo, allTaxId


def extractTaxoFrom_accVSos_ocBDB(acc, allTaxo, BDB):
    # 'acc: os_@#$_oc'
    # allTaxo[acc]['orgName'], allTaxo[acc]['taxId'], allTaxo[acc]['taxoLight'], allTaxo[acc]['DE']
    allTaxo[acc] = {'taxId': '', 'taxoLight': '', 'DE': '', 'orgName': ''}
    os_oc = BDB.get(acc)
    if os_oc:
        os_oc_fld = os_oc.split('_@#$_')
        allTaxo[acc]['orgName'] = os_oc_fld[0]
        allTaxo[acc]['taxoFull'] = os_oc_fld[1]

        return allTaxo[acc]['taxoFull'], allTaxo
    else:
        return '', allTaxo

##############################################################################
#
#            MAIN
#
##############################################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='taxoptimizer.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Parse a blast output report and add NCBI Taxonomy database information in each HSP.")

    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--in", dest="tabfh",
                                 help="Tabulated input file. (Recommended, Blast m8 file)",
                                 metavar="File",
                                 type=file,
                                 default=sys.stdin)
    general_options.add_argument("-o", "--out",
                                 action='store',
                                 dest='outfh',
                                 metavar="File",
                                 type=argparse.FileType('w'),
                                 help='Output file',
                                 default=sys.stdout)
    general_options.add_argument("-b", "--bdb",
                                 dest="bdbfile",
                                 help="",
                                 metavar="File",
                                 required=True
                                 )
    general_options.add_argument("-t", "--bdb_type",
                                 dest="bdbtype",
                                 help="",
                                 type='choice',
                                 default='ncbi',
                                 choices=['ncbi', 'gg', 'silva'],
                                 required=True
                                 )
    general_options.add_argument("-c", "--column",
                                 action='store',
                                 dest='column',
                                 type=int,
                                 help='Column\'s number with ID and database informations for all HSPs',
                                 default=2)
    general_options.add_argument("-s", "--separator",
                                 dest="separator", metavar="str", type=str,
                                 help="Separator in database AC",
                                 default='|')
    general_options.add_argument('-d', '--database', metavar='str',
                                 dest='database',
                                 type=str,
                                 help="Database to use. Supposed that all HSPs match this database.",
                                 )
    general_options.add_argument("-e", "--description",
                                 dest="description",
                                 help="Add database description (DE) in the output",
                                 action='store_true',
                                 default=False,)

    general_options.add_argument("-x", "--splitfile",
                                 dest="splitfile",
                                 help="Only show lines with a taxonomy correspondance in the NCBI taxomomy database. Could be used with -m option.",
                                 action='store_true',
                                 default=False,)
    general_options.add_argument("-f", "--no_taxo_file",
                                 action='store',
                                 dest='notaxofh',
                                 metavar="File",
                                 type=argparse.FileType('w'),
                                 help='Only show lines without a taxonomy correspondance in the NCBI taxomomy database. Could be used with -x option.',)

    golden_options = parser.add_argument_group(title="Golden options", description=None)
    golden_options.add_argument("-m", "--max_cards",
                                action='store',
                                dest='max_cards',
                                type=int,
                                help='Maximum cards number use by Golden2.0 in analyses',
                                default=500)

    args = parser.parse_args()

    cnt_cards = 0
    l_lines = []

    # ===== Tabulated file parsing
    NCBITAXODB_BDB = 'taxodb.bdb'
    GGTAXODB_BDB = 'greengenes_accVosoc.bdb'
    SILVATAXODB_BDB = 'silva_accVosoc.bdb'

    silva_first_pass = True  # open only once silva BDB
    gg_first_pass = True  # open only gg once BDB
    ncbi_osVSocBDB = bdb.DB()

    args.bdbfile = NCBITAXODB_BDB
    try:
        ncbi_osVSocBDB.open(args.bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
    except StandardError, err:
        print >>sys.stderr, TaxOptimizerError("NCBI TaxoDB database open error, %s" % err)
        sys.exit()

    allTaxo = {}
    allTaxId = {}
    try:
        line = args.tabfh.readline()
        lineNb = 1
    except EOFError, err:
        print >>sys.stderr, err
        sys.exit()

    l_cards = ""

    while line:
        description = ''
        fld = line.split()
        if line == '\n':
            line = args.tabfh.readline()
            continue
        try:
            fldcolumn = fld[args.column - 1].split(args.separator)
        except:
            print >>sys.stderr, TaxOptimizerError("Parsing: column error: couldn't parse line: \n%s ...\n --> %s" % (line[0:50], lineNb))
            sys.exit()

        acc = ''
        db = ''
        if len(fldcolumn) == 5:
            # PF8/Pathoquest
            if args.database:
                db = args.database
            else:
                db = fldcolumn[2]
            acc = fldcolumn[3].split('.')[0]
        elif len(fldcolumn) == 2 or len(fldcolumn) == 3:
            # Lionel extraction
            if args.database:
                db = args.database
            else:
                db = fldcolumn[0]
            if len(fldcolumn) == 3 and fldcolumn[1] == '':
                acc = fldcolumn[2].split('.')[0]
            else:
                acc = fldcolumn[1].split('.')[0]
        elif len(fldcolumn) == 1 and args.database:
                db = args.database
                acc = fldcolumn[0]

        if not acc or not db:
            if not args.splitfile:
                print >>args.outfh, line[:-1]
            if args.notaxofh:
                print >>args.notaxofh, line[:-1]
            print >>sys.stderr, TaxOptimizerError("Parsing: acc or db error: %s in line %s" % (fld[args.column - 1], lineNb))

            try:
                line = args.tabfh.readline()
                lineNb += 1
            except EOFError, err:
                print >>sys.stderr, err
                print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))
                sys.exit()
            continue
        elif db not in ['silva', 'gg']:
            l_cards, cnt_cards, l_lines = buildQueryStr(line[:-1], db, acc, l_cards, cnt_cards, l_lines)
            if cnt_cards == args.max_cards:
                allTaxo = doGoldenMulti(allTaxo, l_cards, args.description, allTaxId, ncbi_osVSocBDB)
                printResults(l_lines, allTaxo, args.outfh, args.notaxofh, args.splitfile)
                l_cards = ""
                l_lines = []
                cnt_cards = 0

        else:  # special treatment for silva and gg
            taxonomy = ''
            if acc in allTaxo:
                if 'taxoFull' in allTaxo[acc]:
                    taxonomy = allTaxo[acc]['taxoFull']
                else:
                    taxonomy = allTaxo[acc]['taxoLight']
                description = allTaxo[acc]['DE']
            else:
                taxonomy = ''
                allTaxo[acc] = {'db': db}

                if db == 'silva':
                    if silva_first_pass:
                        silva_first_pass = False
                        silva_accVosocBDB = bdb.DB()
                        silva_accVosocBDBfile = SILVA_TABLE + SILVATAXODB_BDB
                        print >>sys.stderr, silva_accVosocBDBfile
                        try:
                            silva_accVosocBDB.open(silva_accVosocBDBfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
                        except StandardError, err:
                            print >>sys.stderr, TaxOptimizerError("Silva Taxonomy Berkeley open error, %s" % err)
                            sys.exit()
                    taxonomy, allTaxo = extractTaxoFrom_accVSos_ocBDB(acc, allTaxo, silva_accVosocBDB)
                elif db == 'gg':
                    if gg_first_pass:
                        gg_first_pass = False
                        gg_first_pass = False
                        gg_accVosocBDB = bdb.DB()
                        gg_accVosocBDBfile = GG_TABLE + GGTAXODB_BDB

                        try:
                            gg_accVosocBDB.open(gg_accVosocBDBfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
                        except StandardError, err:
                            print >>sys.stderr, TaxOptimizerError("GreenGenes Taxonomy Berkeley database open error, %s" % err)
                            sys.exit()
                    taxonomy, allTaxo = extractTaxoFrom_accVSos_ocBDB(acc, allTaxo, gg_accVosocBDB)

            if taxonomy:
                print >>args.outfh, line[:-1], "\t%s\t%s\t%s" % (allTaxo[acc]['orgName'], taxonomy, allTaxo[acc]['DE'])
            else:
                if args.notaxofh:
                    print >>args.notaxofh, line[:-1]
                if not args.splitFile:
                    print >>args.outfh, line[:-1]

        try:
            line = args.tabfh.readline()
            lineNb += 1
        except EOFError, err:
            print >>sys.stderr, err
            print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))

            sys.exit()
        lineNb += 1

    if cnt_cards != 0:
        allTaxo = doGoldenMulti(allTaxo, l_cards, args.description, allTaxId, ncbi_osVSocBDB)
        printResults(l_lines, allTaxo, args.outfh, args.notaxofh, args.splitfile)

    args.tabfh.close()
    ncbi_osVSocBDB.close()
