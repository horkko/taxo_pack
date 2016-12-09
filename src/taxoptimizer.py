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
from pprint import pprint
from bsddb3 import db as bdb
from bsddb3.db import DBError
from taxodb.errors import Error, GoldenError, TaxOptimizerError, ParserError, ParserWarning
try:
    import Golden
except ImportError as err:
    print(GoldenError("%s\n Install the mandatory program golden (https://github.com/C3BI-pasteur-fr/golden)" %
                      str(err)), file=sys.stderr)
    sys.exit(1)

VERBOSE=False
# What about ['pir', 'pdb', 'tpg', 'tpe', 'tpd', 'prf'] => None
DB_MAPPING = {'sp': 'uniprot', 'sw': 'uniprot', 'tr': 'uniprot', 'trembl': 'uniprot', 'emb': 'embl', 'dbj': 'embl',
              'gb': 'genbank', 'gp': 'genpept', 'ref': 'refseq', 'rdpii': 'rdpii', 'embl_wgs': 'embl_wgs',
               'genbank_wgs': 'genbank_wgs'}



class InputLine:
    # vlegrand@pasteur.fr development
    def __init__(self, orig_line, acc, skip_db):
        self.orig_line = orig_line
        self.acc = acc
        self.skip_db = skip_db


def parseUniprot(input=None, desc=None):
    """
    Parses a UniProt or EMBL like flat file

    :param input: UniProt/EMBL entry
    :type input: str
    :param desc: Description tag value
    type desc: str
    :return: orgName, taxId, taxoLight, description
    :rtype: str, str, str, str
    """
    if not input:
        print(ParserError("No UniProt/EMBL entry given"), file=sys.stderr)
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

    :param input: GenBank entry
    :type input: str
    :param desc: Description tag value
    :type desc: str
    :return: orgName, taxId, taxoLight, description
    :rtype: str, str, str, str
    """
    if not input:
        print(ParserError("No GenBank entry given"), file=sys.stderr)
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

    :param input: Input entry
    :type input: str
    :param desc: Description tag value
    :type desc: str
    :return: orgName, taxId, taxoLight, description
    :rtype: str, str, str, str
    """
    if not input:
        print(ParserError("No input entry given"))
        sys.exit(1)
    if input[:2] == 'ID':
        return parseUniprot(input=input, desc=desc)
    elif input[:5] == 'LOCUS':
        return parseGenbank(input=input, desc=desc)
    return '', '', '', ''

# Builds query input string
#def buildQueryStr(txt_line, db, acc, l_cards, cnt_cards, l_lines):
def buildQueryStr(db=None, acc=None):
    # vlegrand@pasteur.fr development
    if not db or not acc:
        print(TaxOptimizerError("No db or acc given"), file=sys.stderr)
        sys.exit(1)
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
    else:
        skip_db = True
    if skip_db:
        print("[BuildQuery] %s not found" % str(acc))
        return None

    return "%s:%s" % (str(db), str(acc))

##############################################################################
#
#            Golden Multi. Using new version of golden.
#            #vlegrand@pasteur.fr
#
#  l_input : Depending on input_flag, a list of bank:AC bank:locus separated by '\n'
#  so, bank:AC\nbank:locus...
##############################################################################
def doGoldenMulti(taxonomy=None, ids=None, desc=None, taxid=None, bdb=None):
    """
    Performs a golden search

    :param taxonomy: Accession number taxonomy
    :type taxonomy: dict
    :param ids: List of accession/ID to search in database
    :type ids: list
    :param desc: Description type
    :type desc: str
    :param bdb:
    :return:
    """
    if not taxonomy:
        taxonomy = {}
    if not bdb:
        print(GoldenError("No Berkeley db given"), file=sys.stderr)
        sys.exit(1)
    if not len(ids):
        print(GoldenError("No ids to search"), file=sys.stderr)
        sys.exit(1)

    #
    # try:
    #     entry = Golden.access_new(ids)
    #     while entry is not None:
    #         db_acc = lst_input[idx_res]
    #         l_db_acc = db_acc.split(":")
    #         acc = l_db_acc[1]
    #         if acc not in taxonomy:
    #             taxonomy[acc] = {'db': l_db_acc[0]}
    #             taxonomy[acc]['orgName'],\
    #             taxonomy[acc]['taxId'],\
    #             taxonomy[acc]['taxoLight'],\
    #             taxonomy[acc]['DE'] = parse(input=entry, desc=desc)  # orgName, taxId, taxoLight, description
    #             taxonomy, taxid = extractTaxoFrom_osVSocBDB_multi(acc, taxonomy, taxid, bdb)
    #
    #         entry = Golden.access_new(ids)
    #         idx_res += 1
    #     return taxonomy
    # except IOError as err:
    #     print(GoldenError("%s %s" % (str(err), ids)), file=sys.stderr)
    #     sys.exit(1)
    try:
        print("[GOLDEN_MULTI] We have %d entries to search" % len(ids))
        for entry in ids:
            db, acc = entry.split(":")
            if entry == '' or acc in taxonomy:
                continue
            # flat = Golden.access_new(entry)
            flat = Golden.access(db, acc)
            #print("[GOLDEN_MULTI] Results for %s\n%s\n" % (str(entry), str(flat)))
            if flat is not None:
                #print("\tFOUND")
                if acc not in taxonomy:
                    taxonomy[acc] = {'db': db}
                    taxonomy[acc]['orgName'], \
                    taxonomy[acc]['taxId'], \
                    taxonomy[acc]['taxoLight'], \
                    taxonomy[acc]['DE'] = parse(input=flat, desc=desc)
                    taxonomy, taxid = extractTaxoFrom_osVSocBDB_multi(acc, taxonomy, taxid, bdb)
            else:
                pass
                #print("\t** NOT FOUND %s **" % str(flat))
        return taxonomy
    except IOError as err:
        print(GoldenError("%s %s" % (str(err), ids)), file=sys.stderr)
        sys.exit(1)

# display results from taxonomy dictionnary.
# def printResults(lines=None, taxonomy=None, outfile=None, notaxofile=None, splitfile=None):
def printResults(l_lines, taxonomy, outfh, notaxfhout, splitFile):
    if not outfh:
        outfh = sys.stdout

    for li in l_lines:
        print("[%s] ORIGLINE: %s" % (li.acc, str(li.orig_line)), file=outfh)
        taxo = ''
        if not li.skip_db:
            if 'taxoFull' in taxonomy[li.acc]:
                taxo = taxonomy[li.acc]['taxoFull']
            else:
                taxo = taxonomy[li.acc]['taxoLight']
        else:
            if VERBOSE:
                print("SKIP %s" % str(li.orig_line), file=outfh)
        if taxo:
            print("%s\t%s\t%s\t%s" % (str(li.orig_line), str(taxonomy[li.acc]['orgName']), str(taxo),
                                      str(taxonomy[li.acc]['DE'])), file=outfh)
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


def column_analyser(last_field, db):
    """
    Look for accession number in BLAST output match result
    :param last_field: Last field found using separator in tabulated output result file
    :type last_field: list
    :param db: Database we search BLAST hit against (optional)
    :type db: str
    :return: Accession number, Database name
    :rtype: str, str
    """
    acc = ''
    if VERBOSE:
        print("[column_analyzer] last_field=%s (%s)" % (str(last_field), str(db)))
    if len(last_field) == 5:
        if not db:
            db = last_field[2]
        acc = last_field[3].split('.')[0]
    # tr|'D0PRM9|D0PRM9_DUGBV
    elif len(last_field) == 2 or len(last_field) == 3:
        if not db:
            db = last_field[0]
        if len(last_field) == 3 and last_field[1] == '':
            acc = last_field[2].split('.')[0]
        else:
            acc = last_field[1].split('.')[0]
    elif len(last_field) == 1 and db:
        acc = last_field[0]
    # if VERBOSE:
    #     print("[column_analyzer] acc=%s, db=%s" % (acc, db))
    return acc, db


# main_ncbi(tabfh=args.tabfh, outfh=args.outfh, bdb=osVSoc_bdb, column=args.column, sep=args.separator,
#          max_cards=args.max_cards, notaxofh=args.notaxofh, db=args.database, splitfile=args.splitfile,
#          description=args.description)
def main_ncbi(tabfh=None, outfh=None, bdb=None, column=None, separator=None, max_cards=None,
              notaxofh=None, db=None, splitfile=False, description=False):
    if not bdb:
        Error.fatal("No bdb file given")
    if not tabfh:
        Error.fatal("No tabulated (input) file given")
    if not outfh:
        outfh = sys.stdout
    else:
        outfh = open(outfh, 'w')

    taxonomy = {}
    taxids = {}

    try:
        line = tabfh.readline()
        line_number = 1
    except EOFError as err:
        print(TaxOptimizerError("%s" % str(err)), file=sys.stderr)
        sys.exit(1)
    l_cards = ""
    entries = []
    cnt_cards = 0
    l_lines = []
    while line:
        # We can only keep the field of interest instead of keeping all fileds from the complete line
        fields = line.split()
        # Empty line?
        if line == '\n':
            line = tabfh.readline()
            continue
        try:
            id_fields = fields[column - 1].split(separator)
        except Exception as err:
            print(TaxOptimizerError("[%s] Parsing: column error: couldn't parse line: \n%s ...\n --> %s" %
                                    (str(err), line[0:50], str(line_number))), file=sys.stderr)
            outfh.close()
            sys.exit(1)

        acc, db = column_analyser(id_fields, db)

        if not acc or not db:
            if not splitfile:
                print("%s" % str(line[:-1]), file=outfh)
            if notaxofh:
                print("%s" % str(line[:-1]), file=notaxofh)
            print(TaxOptimizerError("Parsing: no acc nor db found in %s using separator=%s (line %s)" %
                                    (fields[column - 1], separator, str(line_number))), file=sys.stderr)
            try:
                line = tabfh.readline()
                line_number += 1
            except EOFError as err:
                print("%s" % str(err), file=sys.stderr)
                print(TaxOptimizerError("at line %s" % str(line_number)), file=sys.stderr)
                outfh.close()
                sys.exit(1)
            continue
        else:
            # l_cards, cnt_cards, l_lines = buildQueryStr(line[:-1], db, acc, l_cards, cnt_cards, l_lines)
            #entry = buildQueryStr(db=db, acc=acc)
            skipdb=True
            if db in DB_MAPPING:
                entry = "%s:%s" % (DB_MAPPING[db], str(acc))
                entries.append(entry)
                cnt_cards += 1
                skipdb=False
                if acc == 'O89815':
                    print("Adding %s" % entry)
            l_lines.append(InputLine(line[:-1], acc, skipdb))
            if cnt_cards == max_cards:
                # taxonomy = doGoldenMulti(taxonomy=taxonomy, ids=l_cards, desc=description, taxid=taxids, bdb=bdb)
                taxonomy = doGoldenMulti(taxonomy=taxonomy, ids=entries, desc=description, taxid=taxids, bdb=bdb)
                # printResults(lines=l_lines, taxonomy=taxonomy, outfile=outfh, notaxofile=notaxofh, spltifile=splitfile)
                printResults(l_lines, taxonomy, outfh, notaxofh, splitfile)
                l_cards = ""
                entries = []
                l_lines = []
                cnt_cards = 0

        try:
            line = tabfh.readline()
            line_number += 1
        except EOFError as err:
            print("%s" % str(err), file=sys.stderr)
            print(TaxOptimizerError("at line %s" % str(line_number)), file=sys.stderr)
            outfh.close()
            sys.exit(1)
        line_number += 1

    if cnt_cards != 0:
        taxonomy = doGoldenMulti(taxonomy=taxonomy, ids=entries, desc=description, taxid=taxids, bdb=bdb)
        if VERBOSE:
            print("** => Golden search done\n%s\n******************" % pprint(str(taxonomy)))
            print("** => l_lines has %d entries" % len(l_lines))
        if 'A6XA53' not in taxonomy:
            print("A6XA53 is not in taxonomy. ABORTED!")
            sys.exit(0)
        # printResults(lines=l_lines, taxonomy=taxonomy, outfile=outfh, notaxofile=notaxofh, splitfile=splitfile)
        printResults(l_lines, taxonomy, outfh, notaxofh, splitfile)
    outfh.close()

def main_gg_silva(tabfh, outfh, accVosocBDB, column, separator,
                  notaxofh=None, db=None, splitfile=False, description=False):

    # allTaxo = {}
    taxonomy = {}
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
            field = fld[column - 1].split(separator)
        except Exception as err:
            print(TaxOptimizerError("[%s] Parsing: column error: couldn't parse line: \n%s ...\n --> %s" %
                                    (str(err), line[0:50], str(lineNb))), file=sys.stderr)
            sys.exit(1)

        acc, db = column_analyser(field, db)
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
            taxo = ''
            if acc in allTaxo:
                if 'taxoFull' in taxonomy[acc]:
                    taxo = taxonomy[acc]['taxoFull']
                else:
                    taxo = taxonomy[acc]['taxoLight']
                if description:
                    DE = allTaxo[acc]['DE']
            else:
                taxo = ''
                allTaxo[acc] = {'db': db}
                taxonomy, allTaxo = extractTaxoFrom_accVSos_ocBDB(acc, allTaxo, accVosocBDB)

            if taxo:
                print("%s\t%s\t%s\t%s" % (line[:-1], taxonomy[acc]['orgName'], taxo, DE), file=outfh)
            else:
                if notaxofh:
                    print("%s", str(line[:-1]), file=notaxofh)
                if not splitfile:
                    print("%s" % str(line[:-1]), file=outfh)

        try:
            line = tabfh.readline()
            lineNb += 1
        except EOFError as err:
            print(TaxOptimizerError("at line %s: %s" % (str(lineNb), str(err))), file=sys.stderr)
            sys.exit(1)
        lineNb += 1


##########################
##         MAIN         ##
##########################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='taxoptimizer.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Parse a blast output report and add NCBI, SILVA or Greengenes"
                                                 "Taxonomy database information in each HSP.")

    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--in",
                                 dest="tabfh",
                                 help="Tabulated input file. (Recommended, Blast m8 file)",
                                 metavar="File name",
                                 type=file,
                                 required=True)
    general_options.add_argument("-o", "--out",
                                 dest='outfh',
                                 type=str,
                                 metavar="File name",
                                 # type=str,
                                 # type=argparse.FileType('w'),
                                 # default=sys.stdout,
                                 default=None,
                                 help='Output file, default prints to STDOUT',
                                 required=False)
    general_options.add_argument("-b", "--bdb",
                                 dest="bdbfile",
                                 help="Berleley DB file generated by taxodb_ncbi"
                                      "(https://github.com/C3BI-pasteur-fr/taxodb_ncbi)"
                                      "or taxo_rrna programs (https://github.com/C3BI-pasteur-fr/taxo_rrna)",
                                 metavar="File",
                                 required=True
                                 )
    general_options.add_argument("-t", "--bdb_type",
                                 dest="bdbtype",
                                 help="Berleley DB file type. Link to -b option.",
                                 type=str,
                                 default='ncbi',
                                 choices=['ncbi', 'greengenes', 'silva'],
                                 required=True
                                 )
    general_options.add_argument("-c", "--column",
                                 action='store',
                                 dest='column',
                                 type=int,
                                 help="Column's number with ID and database information for all HSPs",
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
                                 help="Only show lines with a taxonomy correspondence. Can be used with -m option.",
                                 action='store_true',
                                 default=False,)
    general_options.add_argument("-f", "--no_taxo_file",
                                 dest='notaxofh',
                                 metavar="File name",
                                 type=argparse.FileType('w'),
                                 help='Only show lines without taxonomy correspondence. Can be used with -x option.',
                                 default=False)
    general_options.add_argument('-d', '--database', metavar='str',
                                 dest='database',
                                 type=str,
                                 help="Supposed all blast HSPs match this database."
                                      "ONLY for '-t/--db_type ncbi'.",
                                 default=None,
                                 )
    general_options.add_argument('-v', '--verbose', dest="verbose", action="store_true", default=False,
                                 help="Prints verbose messages")
    golden_options = parser.add_argument_group(title="Golden options", description=None)
    golden_options.add_argument("-m", "--max_cards",
                                action='store',
                                dest='max_cards',
                                type=int,
                                help='Maximum cards number used by Golden >= 2.x in analyses',
                                default=500)

    args = parser.parse_args()
    VERBOSE = args.verbose
    # Only check for GOLDENDATA env variable after all required options are set
    try:
        GOLDENDATA = os.environ['GOLDENDATA']
    except (KeyError, TypeError) as err:
        print(GoldenError('Set the mandatory GOLDENDATA environment variable.\n'
                          'Consult https://github.com/C3BI-pasteur-fr/golden.'),
              file=sys.stderr)
        sys.exit(1)

    # ===== Tabulated file parsing
    if args.bdbtype == 'ncbi':
        try:
            osVSoc_bdb = bdb.DB()
            osVSoc_bdb.open(args.bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
            main_ncbi(tabfh=args.tabfh, outfh=args.outfh, bdb=osVSoc_bdb, column=args.column, separator=args.separator,
                      max_cards=args.max_cards, notaxofh=args.notaxofh, db=args.database, splitfile=args.splitfile,
                      description=args.description)
            # main_ncbi(args.tabfh, args.outfh, osVSoc_bdb, args.column, args.separator, args.max_cards, args.notaxofh,
            #           args.database, args.splitfile, args.description)
            osVSoc_bdb.close()
        except DBError as err:
            print(TaxOptimizerError("Opening %s Berkeley database failed: %s" % (args.bdbtype, str(err))),
                  file=sys.stderr)
            sys.exit(1)

    elif args.bdbtype in ['gg', 'silva']:
        accVosocBDB = bdb.DB()
        try:
            accVosocBDB.open(args.bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
        except StandardError as err:
            print(TaxOptimizerError("Taxonomy Berkeley database open error, %s" % str(err)), file=sys.stderr)
            sys.exit(1)
        # main_gg_silva(args.tabfh, args.outfh, accVosocBDB, args.column, args.separator, args.notaxofh, args.database,
        #               args.splitfile, args.description)
        main_gg_silva(args.tabfh, args.outfh, accVosocBDB, args.column, args.separator, args.notaxofh, args.database,
                      args.splitfile, args.description)
    else:
        print("bdb_type %s not supported" % args.bdbfile, file=sys.stderr)
        sys.exit(1)
    sys.exit(0)
