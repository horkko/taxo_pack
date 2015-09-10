#! /usr/local/bin/python

# Corinne Maufrais
# Institut Pasteur, DSI/CIB
# maufrais@pasteur.fr
#
# golden improvement: vlegrand@pasteur.fr

# version 2.0


import os
import sys
import getopt

import Golden

try:
    GOLDENDATA = os.environ['GOLDENDATA']
except:
    GOLDENDATA = "/local/gensoft2/exe/golden/1.1a/share/golden/db/"
    # GOLDENDATA = "/local/gensoft/share/golden/db"
    os.environ['GOLDENDATA'] = GOLDENDATA

try:
    TABLE = os.environ['TAXOBDBDATA']
except:
    # TABLE = '/local/databases/rel/taxodb/current/bdb/'
    TABLE = '/Users/maufrais/Developpements/TaxoDBNCBI/test/'

try:
    SILVA_TABLE = os.environ['SILVABDBDATA']
except:
    # SILVA_TABLE = '/local/databases/rel/taxodb/current/bdb/'
    SILVA_TABLE = '/Users/maufrais/Developpements/TaxoDB/data/'

try:
    GG_TABLE = os.environ['GGBDBDATA']
except:
    # GG_TABLE = '/local/databases/rel/taxodb/current/bdb/'
    GG_TABLE = '/Users/maufrais/Developpements/TaxoDB/data/'

TAXODB_BDB = 'taxodb.bdb'
NCBITAXODB_BDB = 'taxodb.bdb'
GGTAXODB_BDB = 'greengenes_accVosoc.bdb'
SILVATAXODB_BDB = 'silva_accVosoc.bdb'

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
    return taxonomy,  allTaxo, allTaxId


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

        return allTaxo[acc]['taxoFull'],  allTaxo
    else:
        return '',  allTaxo
##############################################################################
#
#            MAIN
#
##############################################################################


def usage():
    print """
usage: taxoptimizer [options] -i <infile> -o <outfile>

options:

   -i <file> ... Tabulated file.
   -o <file> ... Output file

   -h        ... Print this message and exit.
   -c <int>  ... Column number to parse (default second column: 2)
   -s <car>  ... Separator character (default '|')
   -d <str>  ... Specified Database name for finding taxonomy in only once database.
   -e        ... Add description (DE) in output
   -f <file> ... Extract line without taxonomy from input file
   -x        ... Only write line with taxonomy in output file

   -m <int>  ... Maximum cards number use by Golden2.0 in analyses.
"""


if __name__=='__main__':
    from bsddb3 import db as bdb

    cnt_cards = 0
    l_lines = []

    try:
        TMP_PATH = os.environ['TMPPATH']
    except:
        TMP_PATH = "/tmp"

    #  ===== Command line parser
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "hi:c:s:d:eo:f:xm:", ["help", ])
    except getopt.GetoptError:
        usage()
        sys.exit(0)

    tabFile = None
    outFile = None
    noTaxoFile = None
    splitFile = False
    column = 1
    separator = '|'
    DE = False
    database = None
    max_cards = 500
    for o, v in opts:  # (opt, value)
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("-i", "--in"):
            tabFile = v
        elif o in ("-o", "--out"):
            outFile = v
        elif o in ("-c", "--col"):
            try:
                column = int(v) - 1
            except:
                print >>sys.stderr, TaxOptimizerError("Integer value is mandatory for column number")
                sys.exit(0)
        elif o in ("-d", "--db"):
            database = v
        elif o in ("-s", "--sep"):
            separator = v
        elif o in ("-e", "--desc"):
            DE = True
        elif o in ("-f", "--notaxout"):
            noTaxoFile = v
        elif o in ("-x", "--sep"):
            splitFile = True
        elif o in ("-m", "--max_cards"):
            max_cards = int(v)
        else:
            usage()
            sys.exit(0)

    if not tabFile:
        usage()
        sys.exit(0)

    # more m8Blast_nrprot_light2.txt| ../../src/taxoptimizer.py -o toto -i /dev/stdin

    # ===== Tabulated file parsing
    try:
        tabfhin = open(tabFile)
    except:
        tabfhin = sys.stdin

    try:
        outfh = open(outFile, 'w')
    except:
        outfh = sys.stdout

    notaxfhout = None
    if noTaxoFile:
        notaxfhout = open(noTaxoFile, 'w')

    silva_first_pass = True  # open only once silva BDB
    gg_first_pass = True  # open only gg once BDB
    ncbi_osVSocBDB = bdb.DB()

    ncbi_osVSocBDBfile = TABLE + NCBITAXODB_BDB
    try:
        ncbi_osVSocBDB.open(ncbi_osVSocBDBfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
    except StandardError, err:
        print >>sys.stderr, TaxOptimizerError("NCBI TaxoDB database open error, %s" % err)
        sys.exit()

    allTaxo = {}
    allTaxId = {}
    try:
        line = tabfhin.readline()
        lineNb = 1
    except EOFError, err:
        print >>sys.stderr, err
        sys.exit()

    l_cards = ""

    while line:
        description = ''
        fld = line.split()
        if line == '\n':
            line = tabfhin.readline()
            continue
        try:
            fldInfo = fld[column].split(separator)
        except:
            print >>sys.stderr,  TaxOptimizerError("Parsing: column error: couldn't parse line: \n%s\n --> %s" % (line, lineNb))
            sys.exit()

        acc = ''
        db = ''
        if len(fldInfo) == 5:
            # PF8/Pathoquest
            if database:
                db = database
            else:
                db = fldInfo[2]
            acc = fldInfo[3].split('.')[0]
        elif len(fldInfo) == 2 or len(fldInfo) == 3:
            # Lionel extraction
            if database:
                db = database
            else:
                db = fldInfo[0]
            if len(fldInfo) == 3 and fldInfo[1] == '':
                acc = fldInfo[2].split('.')[0]
            else:
                acc = fldInfo[1].split('.')[0]
        elif len(fldInfo) == 1 and database:
                db = database
                acc = fldInfo[0]

        if not acc or not db:
            if not splitFile:
                print >>outfh, line[:-1]
            if noTaxoFile:
                print >>notaxfhout,  line[:-1]
            print >>sys.stderr, TaxOptimizerError("Parsing: acc or db error: :%s in line %s" % (fld[column], lineNb))

            try:
                line = tabfhin.readline()
                lineNb += 1
            except EOFError, err:
                print >>sys.stderr, err
                print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))
                sys.exit()
            continue
        elif db not in ['silva', 'gg']:
            l_cards, cnt_cards, l_lines = buildQueryStr(line[:-1], db, acc, l_cards, cnt_cards, l_lines)
            if cnt_cards == max_cards:
                allTaxo = doGoldenMulti(allTaxo, l_cards, DE, allTaxId, ncbi_osVSocBDB)
                printResults(l_lines, allTaxo, outfh, notaxfhout, splitFile)
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
                print >>outfh, line[:-1], "\t%s\t%s\t%s" % (allTaxo[acc]['orgName'], taxonomy, allTaxo[acc]['DE'])
            else:
                if noTaxoFile:
                    print >>notaxfhout, line[:-1]
                if not splitFile:
                    print >>outfh, line[:-1]

        try:
            line = tabfhin.readline()
            lineNb += 1
        except EOFError, err:
            print >>sys.stderr, err
            print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))

            sys.exit()
        lineNb += 1

    if cnt_cards != 0:
        allTaxo = doGoldenMulti(allTaxo, l_cards, DE, allTaxId, ncbi_osVSocBDB)
        printResults(l_lines, allTaxo, outfh, notaxfhout, splitFile)

    tabfhin.close()
    ncbi_osVSocBDB.close()
