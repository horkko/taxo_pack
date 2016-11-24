#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# corinne.maufrais@pasteur.fr
# Updated by Emmanuel Quevillon, tuco@pasteur.fr
# version 2.1
from __future__ import print_function
import os
import sys
import argparse
from bsddb3 import db as bdb
from taxodb.errors import GoldenError, TaxOptimizerError, ParserError, ParserWarning
try:
    import Golden
except ImportError as err:
    print(GoldenError("%s\n Install the mandatory program golden (https://github.com/C3BI-pasteur-fr/golden)" %
                      str(err)), file=sys.stderr)
    sys.exit(1)

try:
    GOLDENDATA = os.environ['GOLDENDATA']
except KeyError as err:
    print(GoldenError('Set the mandatory GOLDENDATA environment variable.\n'
                      'Consult https://github.com/C3BI-pasteur-fr/golden.' % str(err),
                      file=sys.stderr))
    sys.exit(1)

class InputLine:
    # vlegrand@pasteur.fr development
    def __init__(self, orig_line, acc, skip_db):
        self.orig_line = orig_line
        self.acc = acc
        self.skip_db = skip_db


def parseUniprot(input=None, desc=None):
    """
    Parses a UniProt or EMBL like flat file

    :param input: Input file to parse
    :type input: str
    :param desc: Description tag value
    type desc: str
    :return: orgName, taxId, taxoLight, description
    :rtype: str, str, str, str
    """
    if not input:
        print(ParserError("No UniProt/EMBL file given"), file=sys.stderr)
        sys.exit(1)
    description = ''
    taxId = ''
    taxoLight = ''
    orgName = ''
    lineFld = input.split('\n')
    vuOS = vuOC = vuOCXX = False
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:2].strip()
        value = line[5:]
        value = value.replace('\n', '')
        if desc and tag == desc:
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


def parseGenbank(input=None, desc=None):
    """
    Parses GenBank like flat file

    :param input: GenBank input file
    :type input: str
    :param desc: Description tag value
    :type desc: str
    :return: orgName, taxId, taxoLight, description
    :rtype: str, str, str, str
    """
    if not input:
        print(ParserError("No GenBank file given"), file=sys.stderr)
        sys.exit(1)
    description = ''
    taxId = ''
    taxoLight = ''
    orgName = ''
    taxoL = False
    lineFld = input.split('\n')
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:12].strip()
        value = line[12:]
        value = value.replace('\n', '')
        if desc and tag == 'DEFINITION':
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


def parse(input=None, desc=None):
    """
    Parses db flat file (UniProt, EMBL, GenBank) into a DBRecord object

    :param input: Input file
    :type input: str
    :param desc: Description tag value
    :type desc: str
    :return: orgName, taxId, taxoLight, description
    :rtype: str, str, str, str
    """
    if not input:
        print(ParserError("No input file gven"))
        sys.exit(1)
    if input[:2] == 'ID':
        return parseUniprot(input=input, desc=desc)
    elif input[:5] == 'LOCUS':
        return parseGenbank(input=input, desc=desc)
    return '', '', '', ''

##############################################################################
#
#            Golden: Keep this for compatibility and performance testing.
#
##############################################################################


def doGolden(db=None, acc=None, desc=None):
    """
    Performs a golden search

    :param db: Database to search
    :type db: str
    :param acc: Accession/ID to search in database
    :type acc: str
    :param desc: Description type
    :type desc: str
    :return:
    """
    # ########################## db ref
    if not db:
        print(GoldenError("No database to search given"))
        sys.exit(1)
    if not acc:
        print(GoldenError("No Acc/ID to search given"))
        sys.exit(1)

    if db in ['sp', 'sw', 'swissprot', 'tr', 'trembl']:
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
        flatFile = Golden.access(db, acc)
    except IOError as err:
        print("%s, %s, %s" % (str(err), db, acc), file=sys.stderr)
        sys.exit(1)
    if flatFile:
        orgName, taxId, taxoLight, description = parse(input=flatFile, desc=desc)  # orgName, taxId, taxoLight, description
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


def doGoldenMulti(allTaxo, l_cards, DE, allTaxId, osVSoc_bdb):
    idx_res = 0
    lst_input = l_cards.split("\n")
    try:
        flatFile = Golden.access_new(l_cards)
        while (flatFile is not None):
            db_acc = lst_input[idx_res]
            l_db_acc = db_acc.split(":")
            acc = l_db_acc[1]
            if acc not in allTaxo:
                allTaxo[acc] = {'db': l_db_acc[0]}
                allTaxo[acc]['orgName'], allTaxo[acc]['taxId'], allTaxo[acc]['taxoLight'], allTaxo[acc]['DE'] = parse(flatFile, DE)  # orgName, taxId, taxoLight, description
                allTaxo, allTaxId = extractTaxoFrom_osVSocBDB_multi(acc, allTaxo, allTaxId, osVSoc_bdb)

            flatFile = Golden.access_new(l_cards)
            idx_res += 1
        return allTaxo
    except IOError as err:
        print("%s %s" % (str(err), l_cards), file=sys.stderr)
        sys.exit(1)

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
            print("%s\t%s\t%s\t%s" % (str(li.orig_line), str(allTaxo[li.acc]['orgName']), str(taxonomy),
                                      str(allTaxo[li.acc]['DE'])), file=outfh)
        else:
            if notaxfhout:
                print("%s" % str(li.orig_line), file=notaxfhout)
            if not splitFile:
                print("%s" % str(li.orig_line), file=outfh)

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


def column_analyser(fldcolumn, db):
    acc = ''
    if len(fldcolumn) == 5:
        if not db:
            db = fldcolumn[2]
        acc = fldcolumn[3].split('.')[0]
    elif len(fldcolumn) == 2 or len(fldcolumn) == 3:
        if not db:
            db = fldcolumn[0]
        if len(fldcolumn) == 3 and fldcolumn[1] == '':
            acc = fldcolumn[2].split('.')[0]
        else:
            acc = fldcolumn[1].split('.')[0]
    elif len(fldcolumn) == 1 and db:
        acc = fldcolumn[0]
    return acc, db


def main_ncbi(tabfh, outfh, osVSoc_bdb, column, separator, max_cards,
              notaxofh=None, db=None, splitfile=False, description=False):
    allTaxo = {}
    allTaxId = {}
    try:
        line = tabfh.readline()
        lineNb = 1
    except EOFError as err:
        print("%s" % str(err), file=sys.stderr)
        sys.exit(1)
    l_cards = ""
    cnt_cards = 0
    l_lines = []
    while line:
        fld = line.split()
        if line == '\n':
            line = tabfh.readline()
            continue
        try:
            fldcolumn = fld[column - 1].split(separator)
        except Exception as err:
            print(TaxOptimizerError("[%s] Parsing: column error: couldn't parse line: \n%s ...\n --> %s" %
                                    (str(err), line[0:50], str(lineNb))), file=sys.stderr)
            sys.exit(1)

        acc, db = column_analyser(fldcolumn, db)

        if not acc or not db:
            if not splitfile:
                print("%s" % str(line[:-1]), file=outfh)
            if notaxofh:
                print("%s" % str(line[:-1]), file=notaxofh)
            print(TaxOptimizerError("Parsing: no acc and db in %s with separator=%s (line %s)" %
                                    (fld[column - 1], separator, str(lineNb))), file=sys.stderr)
            try:
                line = tabfh.readline()
                lineNb += 1
            except EOFError as err:
                print("%s" % str(err), file=sys.stderr)
                print(TaxOptimizerError("at line %s" % str(lineNb)), file=sys.stderr)
                sys.exit(1)
            continue
        elif db not in ['silva', 'gg']:
            l_cards, cnt_cards, l_lines = buildQueryStr(line[:-1], db, acc, l_cards, cnt_cards, l_lines)
            if cnt_cards == max_cards:
                allTaxo = doGoldenMulti(allTaxo, l_cards, description, allTaxId, osVSoc_bdb)
                printResults(l_lines, allTaxo, outfh, notaxofh, splitfile)
                l_cards = ""
                l_lines = []
                cnt_cards = 0

        try:
            line = tabfh.readline()
            lineNb += 1
        except EOFError as err:
            print("%s" % str(err), file=sys.stderr)
            print(TaxOptimizerError("at line %s" % str(lineNb)), file=sys.stderr)
            sys.exit(1)
        lineNb += 1

    if cnt_cards != 0:
        allTaxo = doGoldenMulti(allTaxo, l_cards, description, allTaxId, osVSoc_bdb)
        printResults(l_lines, allTaxo, outfh, notaxofh, splitfile)


def main_gg_silva(tabfh, outfh, accVosocBDB, column, separator,
                  notaxofh=None, db=None, splitfile=False, description=False):
    allTaxo = {}
    try:
        line = tabfh.readline()
        lineNb = 1
    except EOFError as err:
        print("%s" % str(err), file=sys.stderr)
        sys.exit(1)
    DE = ''
    while line:
        fld = line.split()
        if line == '\n':
            line = tabfh.readline()
            continue
        try:
            fldcolumn = fld[column - 1].split(separator)
        except Exception as err:
            print(TaxOptimizerError("[%s] Parsing: column error: couldn't parse line: \n%s ...\n --> %s" %
                                    (str(err), line[0:50], str(lineNb))), file=sys.stderr)
            sys.exit(1)

        acc, db = column_analyser(fldcolumn, db)
        if not acc or not db:
            if not splitfile:
                print("%s" % str(line[:-1]), file=outfh)
            if notaxofh:
                print("%s" % str(line[:-1]), file=notaxofh)
            print(TaxOptimizerError("Parsing: no acc and db in %s with separator=%s (line %s)" %
                                    (fld[column - 1], separator, str(lineNb))), file=sys.stderr)
            try:
                line = tabfh.readline()
                lineNb += 1
            except EOFError as err:
                print(TaxOptimizerError("at line %s" % (lineNb)), file=sys.stderr)
                sys.exit(1)
            continue
        else:
            taxonomy = ''
            if acc in allTaxo:
                if 'taxoFull' in allTaxo[acc]:
                    taxonomy = allTaxo[acc]['taxoFull']
                else:
                    taxonomy = allTaxo[acc]['taxoLight']
                if description:
                    DE = allTaxo[acc]['DE']
            else:
                taxonomy = ''
                allTaxo[acc] = {'db': db}
                taxonomy, allTaxo = extractTaxoFrom_accVSos_ocBDB(acc, allTaxo, accVosocBDB)

            if taxonomy:
                print("%s\t%s\t%s\t%s" % (line[:-1], allTaxo[acc]['orgName'], taxonomy, DE), file=outfh)
            else:
                if notaxofh:
                    print("%s", str(line[:-1]), file=notaxofh)
                if not splitfile:
                    print("%s" % str(line[:-1]), file=outfh)

        try:
            line = tabfh.readline()
            lineNb += 1
        except EOFError as err:
            print("%s" % str(err), file=sys.stderr)
            print(TaxOptimizerError("at line %s" % str(lineNb)), file=sys.stderr)
            sys.exit(1)
        lineNb += 1


##############################################################################
#
#            MAIN
#
##############################################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='taxoptimizer.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Parse a blast output report and add NCBI, SILVA or Greengenes Taxonomy database information in each HSP.")

    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--in", dest="tabfh",
                                 help="Tabulated input file. (Recommended, Blast m8 file)",
                                 metavar="File",
                                 type=file,
                                 required=True)
    general_options.add_argument("-o", "--out",
                                 action='store',
                                 dest='outfh',
                                 metavar="File",
                                 type=argparse.FileType('w'),
                                 help='Output file',
                                 required=True)
    general_options.add_argument("-b", "--bdb",
                                 dest="bdbfile",
                                 help="Berleley DB file generated by taxodb_ncbi (https://github.com/C3BI-pasteur-fr/taxodb_ncbi) or taxo_rrna programs (https://github.com/C3BI-pasteur-fr/taxo_rrna)",
                                 metavar="File",
                                 required=True
                                 )
    general_options.add_argument("-t", "--bdb_type",
                                 dest="bdbtype",
                                 help="Berleley DB file type. Link to -b option.",
                                 type=str,
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

    general_options.add_argument("-e", "--description",
                                 dest="description",
                                 help="Add database description (DE) in the output.",
                                 action='store_true',
                                 default=False,)
    general_options.add_argument("-x", "--splitfile",
                                 dest="splitfile",
                                 help="Only show lines with a taxonomy correspondance. Could be used with -m option.",
                                 action='store_true',
                                 default=False,)
    general_options.add_argument("-f", "--no_taxo_file",
                                 action='store',
                                 dest='notaxofh',
                                 metavar="File",
                                 type=argparse.FileType('w'),
                                 help='Only show lines without a taxonomy correspondance. Could be used with -x option.',)
    general_options.add_argument('-d', '--database', metavar='str',
                                 dest='database',
                                 type=str,
                                 help="Supposed that all blast HSPs match this database. Not used for Silva and Greengenes databases",
                                 default=None,
                                 )
    golden_options = parser.add_argument_group(title="Golden options", description=None)
    golden_options.add_argument("-m", "--max_cards",
                                action='store',
                                dest='max_cards',
                                type=int,
                                help='Maximum cards number used by Golden2.0 in analyses',
                                default=500)

    args = parser.parse_args()

    # ===== Tabulated file parsing
    if args.bdbtype == 'ncbi':
        osVSoc_bdb = bdb.DB()
        try:
            osVSoc_bdb.open(args.bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
        except StandardError as err:
            print(TaxOptimizerError("NCBI TaxoDB database open error, %s" % str(err)), file=sys.stderr)
            sys.exit(1)

        main_ncbi(args.tabfh, args.outfh, osVSoc_bdb, args.column, args.separator, args.max_cards, args.notaxofh, args.database, args.splitfile, args.description)
        osVSoc_bdb.close()
    elif args.bdbtype in ['gg', 'silva']:
        accVosocBDB = bdb.DB()
        try:
            accVosocBDB.open(args.bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
        except StandardError as err:
            print(TaxOptimizerError("Taxonomy Berkeley database open error, %s" % str(err)), file=sys.stderr)
            sys.exit(1)
        main_gg_silva(args.tabfh, args.outfh, accVosocBDB, args.column, args.separator, args.notaxofh, args.database, args.splitfile, args.description)
    sys.exit(0)