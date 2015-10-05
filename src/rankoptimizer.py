#! /local/gensoft2/adm/bin/python

# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# maufrais@pasteur.fr
#

# version 1.6

import os
import sys
import argparse


try:
    LIB = os.environ['RANKOPTIMIZERLIB']
except:
    LIB = '/usr/local/bin'
    LIB = '/Users/maufrais/Developpements2/taxo_pack/lib'

if LIB not in sys.path:
    sys.path.append(LIB)

import rankoptimizerlib


class RankOptimizerError:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[rankoptimizer] " + self.err


def clean_taxo(sbjct_taxonomy, clean_cellular):
    # -r        ... Clean genbank taxonomy error format (default:false). Very very slow.
    fld1 = sbjct_taxonomy.split('Viruses')
    if len(fld1) > 1 and fld1[0] != "''":
        if fld1[0] != "cellular organisms; " or (clean_cellular and fld1[0] == "cellular organisms; "):
            sbjct_taxonomy = 'Viruses' + ''.join(fld1[1:])

    fld2 = sbjct_taxonomy.split('Bacteria')
    if len(fld2) > 1 and fld2[0] != "''":
        if fld2[0] != "cellular organisms; " or (clean_cellular and fld2[0] == "cellular organisms; "):
            sbjct_taxonomy = 'Bacteria' + ''.join(fld2[1:])

    fld3 = sbjct_taxonomy.split('Eukaryota')
    if len(fld3) > 1 and fld3[0] != "''":
        if fld3[0] != "cellular organisms; " or (clean_cellular and fld3[0] == "cellular organisms; "):
            sbjct_taxonomy = 'Eukaryota' + ''.join(fld3[1:])

    fld4 = sbjct_taxonomy.split('Archaea')
    if len(fld4) > 1 and fld4[0] != "''":
        if fld4[0] != "cellular organisms; " or (clean_cellular and fld4[0] == "cellular organisms; "):
            sbjct_taxonomy = 'Archaea' + ''.join(fld4[1:])

    fld5 = sbjct_taxonomy.split('other sequences')
    if len(fld5) > 1 and fld5[0] != "''":
        if fld5[0] != "cellular organisms; " or (clean_cellular and fld5[0] == "cellular organisms; "):
            sbjct_taxonomy = 'other sequences' + ''.join(fld5[1:])

    fld6 = sbjct_taxonomy.split('unclassified sequences')
    if len(fld6) > 1 and fld6[0] != "''":
        if fld6[0] != "cellular organisms; " or (clean_cellular and fld6[0] == "cellular organisms; "):
            sbjct_taxonomy = 'unclassified sequences' + ''.join(fld6[1:])

    return sbjct_taxonomy


def extract_tot_rank(taxonomy):
    rank_list = ['superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass',
                 'infraclass', 'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily', 'tribe',
                 'subtribe', 'genus', 'subgenus', 'species_group', 'species_subgroup', 'species', 'subspecies', 'varietas', 'forma']
    fld = taxonomy.split(';')
    all_nodes = []
    for node in fld:
        rank = ''
        fld_node = node.strip().split('(')
        name = fld_node[0].strip()
        if len(fld_node) == 2:
            rank = fld_node[1].split(')')[0].strip()
            rank = rank.replace(' ', '_')
        elif len(fld_node) > 2:
            rank = fld_node[-1].split(')')[0].strip()
            name = fld_node[0] + fld_node[1].split(')')[0].strip()
            rank = rank.replace(' ', '_')
        if rank and rank not in rank_list:
            name += '_' + rank
            rank = ''
        if name:
            all_nodes.append((name, rank))

    return all_nodes


def insert_taxo_tot(taxo_tree, query, pos_line, taxonomy, rank_in_tree, identical_read):
    all_nodes = extract_tot_rank(taxonomy)
    if taxo_tree.name == 'root':
        if identical_read and '_' in query:  # # supposed name_n
            n = query.split('_')[-1]
            try:
                nb = int(n)
            except:
                nb = 1
            taxo_tree.nb_querys += nb
        else:
            taxo_tree.nb_querys += 1

    try:
        name, rank = all_nodes[0]
    except:
        return taxo_tree

    if not taxo_tree.has_child(name):
        if rank_in_tree:
            taxo_tree.add_child(name, rank)
        else:
            taxo_tree.add_child(name, '')
    child = taxo_tree.get_child(name)

    # child.add_queries((query,  pos_line) )

    if identical_read and '_' in query:  # # supposed name_n
        n = query.split('_')[-1]
        try:
            nb = int(n)
        except:
            nb = 1
        child.nb_querys += nb
    else:
        child.nb_querys += 1
    for name, rank in all_nodes[1:-1]:
        parent = child
        if name:
            if not parent.has_child(name):
                if rank_in_tree:
                    parent.add_child(name, rank)
                else:
                    parent.add_child(name, '')
            child = parent.get_child(name)
            # child.add_queries((query,  pos_line))
            if identical_read and '_' in query:  # # supposed name_n
                n = query.split('_')[-1]
                try:
                    nb = int(n)
                except:
                    nb = 1
                child.nb_querys += nb
            else:
                child.nb_querys += 1
            if not child.has_rank() and rank_in_tree and rank:
                child.rank = rank

    # add query name, only on the leaf
    name, rank = all_nodes[-1]
    parent = child
    if name:
        if not parent.has_child(name):
            if rank_in_tree:
                parent.add_child(name, rank)
            else:
                parent.add_child(name, '')
        child = parent.get_child(name)
        child.add_queries((query,  pos_line))
        if identical_read and '_' in query:  # # supposed name_n
            n = query.split('_')[-1]
            try:
                nb = int(n)
            except:
                nb = 1
            child.nb_querys += nb
        else:
            child.nb_querys += 1
        if not child.has_rank() and rank_in_tree and rank:
            child.rank = rank
    return taxo_tree


def insert_taxo_tot_delta(taxo_tree, query, pos_line, taxonomy, rank_in_tree, identical_read):  # ## idem que insert_taxo_tot why ???????
    all_nodes = extract_tot_rank(taxonomy)
    if taxo_tree.name == 'root':
        if identical_read and '_' in query:  # # supposed name_n
            n = query.split('_')[-1]
            try:
                nb = int(n)
            except:
                nb = 1
            taxo_tree.nb_querys += nb
        else:
            taxo_tree.nb_querys += 1

    try:
        name, rank = all_nodes[0]
    except:
        return taxo_tree

    if not taxo_tree.has_child(name):
        if rank_in_tree:
            taxo_tree.add_child(name, rank)
        else:
            taxo_tree.add_child(name, '')
    child = taxo_tree.get_child(name)

    # child.add_queries((query,  pos_line) )
    if identical_read and '_' in query:  # # supposed name_n
        n = query.split('_')[-1]
        try:
            nb = int(n)
        except:
            nb = 1
        child.nb_querys += nb
    else:
        child.nb_querys += 1
    for name, rank in all_nodes[1:-1]:
        parent = child
        if name:
            if not parent.has_child(name):
                if rank_in_tree:
                    parent.add_child(name, rank)
                else:
                    parent.add_child(name, '')
            child = parent.get_child(name)
            # child.add_queries((query,  pos_line))
            if identical_read and '_' in query:  # # supposed name_n
                n = query.split('_')[-1]
                try:
                    nb = int(n)
                except:
                    nb = 1
                child.nb_querys += nb
            else:
                child.nb_querys += 1
            if not child.has_rank() and rank_in_tree and rank:
                child.rank = rank

    # add query name, only on the leaf
    name, rank = all_nodes[-1]
    parent = child
    if name:
        if not parent.has_child(name):
            if rank_in_tree:
                parent.add_child(name, rank)
            else:
                parent.add_child(name, '')
        child = parent.get_child(name)
        child.add_queries((query, pos_line))
        if identical_read and '_' in query:  # # supposed name_n
            n = query.split('_')[-1]
            try:
                nb = int(n)
            except:
                nb = 1
            child.nb_querys += nb
        else:
            child.nb_querys += 1
        if not child.has_rank() and rank_in_tree and rank:
            child.rank = rank

    return taxo_tree


def reduce_by_delta(one_query, query, delta_score):

    # query_infos[query]['delta'] = [(hsp_score, sbjct_taxonomy, pos_line)]
    # query_infos[query]['max_score'] = (hsp_score, sbjct_taxonomy, pos_line)

    # kept only bline in a delta range
    indx = []
    for i in range(len(one_query['delta'])):
        if one_query['delta'][i][0] > (one_query['max_score'][0] * (1 + (delta_score*0.01))):
            indx.append(i)
    n = 0
    for j in indx:
        j -= n
        one_query['delta'].pop(j)
        n += 1

    return one_query


def lca(tree):
    if tree.hasOneChild():
        tree = tree.childs[0]
        return lca(tree)
    else:
        return tree.name

##############################################################################
#
#            MAIN
#
##############################################################################


def usage():
    print """
    rankoptimizer use Krona 2.1 an interactive metagenomic visualization tool in a Web browser.  (http://sourceforge.net/p/krona/home/krona/):
    Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.

"""

if __name__ == '__main__':
    # ===== Command line parser
    parser = argparse.ArgumentParser(prog='rankoptimizer.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Statistical analyses of the NCBI Taxonomy informations contained in a taxoptimizer file.\
                                     Output file could be visualized by Krona 2.1 an interactive metagenomic visualization in a \
                                     Web browser.")

    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--in", dest="tabfh",
                                 help="Tabulated input file. Blast report with NCBI Taxonomy database informations from taxoptimizer programm.",
                                 metavar="File",
                                 required=True)

    general_options.add_argument("-c", "--taxcolumn",
                                 action='store',
                                 dest='taxcolumn',
                                 type=int,
                                 help='First column with NCBI Taxonomy informations.',
                                 default=14)
    general_options.add_argument("-C", "--scorecolumn",
                                 action='store',
                                 dest='scorecolumn',
                                 type=int,
                                 help='Column\'s number with HSP scores',
                                 default=12)

    output_options = parser.add_argument_group(title="Output options", description=None)
    output_options.add_argument("-k", "--krona",
                                action='store',
                                dest='kronafh',
                                metavar="File",
                                type=argparse.FileType('w'),
                                help='xml output file with Krona 2.1 Specification.',)
    output_options.add_argument("-t", "--text",
                                action='store',
                                dest='textfh',
                                metavar="File",
                                type=argparse.FileType('w'),
                                help='taxonomy abundance with dendrogram tree representation',)
    output_options.add_argument("-v", "--htmlx",
                                action='store',
                                dest='htmlxmlfh',
                                metavar="File",
                                type=argparse.FileType('w'),
                                help='html output with Krona 2.0 specification and Krona 2.1 javascript library. xml style',)
    output_options.add_argument("-V", "--htmlj",
                                action='store',
                                dest='htmljsonfh',
                                metavar="File",
                                type=argparse.FileType('w'),
                                help='html output with Krona 2.0 specification and Krona 2.1 javascript library. json style',)
    output_options.add_argument("-j", "--json",
                                action='store',
                                dest='jsonfh',
                                metavar="File",
                                type=argparse.FileType('w'),
                                help='json output with Krona 2.1 specification',)
    output_options.add_argument("-p", "--dump",
                                action='store',
                                dest='dumpfh',
                                metavar="File",
                                type=argparse.FileType('w'),
                                help='Python dump output with Krona 2.1 specification.',)
    output_options.add_argument("-a", "--lca",
                                dest="lca",
                                help="Report lowest common ancestor of the taxonomic abundance",
                                action='store_true',
                                default=False,)

    specific_options = parser.add_argument_group(title="Specific options", description=None)

    specific_options.add_argument("-u", "--url",
                                  dest="kronahome", metavar="str", type=str,
                                  help="KronaTools-2.1 http server address",
                                  default='http://krona.sourceforge.net')

    specific_options.add_argument('-s', '--local',
                                  dest='jslocal',
                                  action='store_true',
                                  default=False,
                                  help="Local Krona javascript.",
                                  )
    specific_options.add_argument("-d", "--delta",
                                  dest="delta",
                                  help="Accepted percentage variation of HSP score value.",
                                  default=0,
                                  type=int,)
    specific_options.add_argument("-l", "--clean",
                                  dest="clean",
                                  help="Remove 'cellular organism' from taxonomy.",
                                  action='store_true',
                                  default=False,)
    specific_options.add_argument("-R", "--rank",
                                  dest="rank",
                                  help="""Report rank information in xml krona (-k) and html krona (-v) formats.
                                  By default, "superkingdom, phylum, ..." information are not analyzed.
                                  Note: If two queries have the same taxonomy but only one has rank information, queries are in differents nodes into the tree.""",
                                  action='store_true',
                                  default=False,)
    specific_options.add_argument("-U", "--identical",
                                  dest="identical",
                                  help="""Reintroduce abundance of identical reads in the analyze.
                                  Identical reads are treat as unique sequence. The number of identical reads is put at the end of the query name (separated by '_').""",
                                  action='store_true',
                                  default=False,)

    args = parser.parse_args()

    # ### args.kronahome = 'http://bioweb2.pasteur.fr/krona/src' ==> Ne fonctionne plus

    query_infos = {}
    taxcolumn = args.taxcolumn - 1
    scorecolumn = args.scorecolumn - 1
    try:
        pos_line = args.tabfh.tell()
        line = args.tabfh.readline()
        blast_line = 1
    except IOError, err:
        print >>sys.stderr, RankOptimizerError("readline error: %s line:%s" % (err, blast_line))
        sys.exit()

    while line:
        if line == '\n':
            try:
                pos_line = args.tabfh.tell()
                line = args.tabfh.readline()
                blast_line += 1
            except IOError, err:
                print >>sys.stderr, RankOptimizerError("readline error:%s line:%s" % (err, blast_line))
                sys.exit()
            continue
        fld = line.split('\t')
        try:
            query = fld[0]
        except StandardError, err:
            print >>sys.stderr, RankOptimizerError("query error:%s line:%s" % (err, blast_line))
            sys.exit()

        try:
            if fld[scorecolumn] and fld[scorecolumn][0] == 'e':
                num = '1' + fld[scorecolumn]
            else:
                num = fld[scorecolumn]
            hsp_score = float(num)
        except StandardError, err:
            print >>sys.stderr, RankOptimizerError("score column number error:-C %s in line %s:%s" % (scorecolumn+1, blast_line, err))
            sys.exit()
        try:
            sbjct_taxonomy = fld[taxcolumn]
            sbjct_taxonomy = True
        except:
            sbjct_taxonomy = False

        # {'hsp_score':hsp_score, 'sbjct_taxonomy':sbjct_taxonomy, 'pos_line':pos_line }
        # ===> (hsp_score,  sbjct_taxonomy, pos_line)
        # ##### Analyze and stock query
        if query in query_infos:
            if not query_infos[query]['max_score'][1] or (sbjct_taxonomy and hsp_score > query_infos[query]['max_score'][0]):
                query_infos[query]['max_score'] = (hsp_score, sbjct_taxonomy, pos_line)
            elif (args.delta and sbjct_taxonomy) and (hsp_score > (query_infos[query]['max_score'][0] * (1 - (args.delta*0.01)))):
                # stock all at one moment. reduce after.
                if 'delta' in query_infos[query]:
                    query_infos[query]['delta'].append((hsp_score, sbjct_taxonomy, pos_line))
                else:
                    query_infos[query]['delta'] = [(hsp_score, sbjct_taxonomy, pos_line)]
            # else:
            #    print ''
        else:
            query_infos[query] = {'max_score': (hsp_score, sbjct_taxonomy, pos_line)}

        try:
            pos_line = args.tabfh.tell()
            line = args.tabfh.readline()
            blast_line += 1
        except IOError, err:
            print >>sys.stderr, RankOptimizerError("%s line:%s" % (err, blast_line))
            sys.exit()

    # tabfhin.close()
    if args.delta:
        for query in query_infos.keys():
            if 'delta' in query_infos[query] and len(query_infos[query]['delta']) > 1:
                # avant= len(query_infos[query]['delta'])
                query_infos[query] = reduce_by_delta(query_infos[query], query, args.delta)

    # ## Construct tree structure
    taxo_tree = rankoptimizerlib.Taxon('root')
    # print >>sys.stderr, 'beginning insert taxo', time.strftime("%y/%m/%d %H:%M:%S" , time.localtime(time.time()))

    # p = psutil.Process(os.getpid())
    for query, infos in query_infos.items():
        pos_line = infos['max_score'][2]
        args.tabfh.seek(pos_line)
        try:
            line = args.tabfh.readline()
            fld = line.split('\t')
        except IOError, err:
            print >>sys.stderr, RankOptimizerError("%s line:%s" % (err, line))
            sys.exit()
        try:
            sbjct_taxonomy = fld[taxcolumn].strip()  # ## taxonomy begin at the 13 column separated by \t
            # ###### Clean taxonomy
            if sbjct_taxonomy[-1] == '.':
                sbjct_taxonomy = sbjct_taxonomy[:-1] + ';'
            if args.clean:
                sbjct_taxonomy = sbjct_taxonomy.replace('cellular organisms ;', '').strip()
                sbjct_taxonomy = sbjct_taxonomy.replace('cellular organisms;', '').strip()
            # if repairGB and sbjct_taxonomy:### slow
            #    sbjct_taxonomy = clean_taxo(sbjct_taxonomy)

        except:
            print >>sys.stderr, 'no taxo found for %s' % query
            continue

        taxo_tree = insert_taxo_tot(taxo_tree, query, pos_line, sbjct_taxonomy, args.rank, args.identical)

        if 'delta' in infos:
            all_delta_taxo = [sbjct_taxonomy]
            for qd in infos['delta']:
                pos_line = qd[2]
                args.tabfh.seek(pos_line)
                blast_line = args.tabfh.readline()
                fld = blast_line.split('\t')
                try:
                    sbjct_taxonomy = fld[taxcolumn].strip()
                    # ###### Clean taxonomy
                    if sbjct_taxonomy[-1] == '.':
                        sbjct_taxonomy = sbjct_taxonomy[:-1] + ';'
                    if args.clean:
                        sbjct_taxonomy = sbjct_taxonomy.replace('cellular organisms ;', '').strip()
                        sbjct_taxonomy = sbjct_taxonomy.replace('cellular organisms;', '').strip()
                    # if repairGB and sbjct_taxonomy:### slow
                    #    sbjct_taxonomy = clean_taxo(sbjct_taxonomy)
                except:
                    print >>sys.stderr, 'no taxo found for %s' % query
                    continue
                if sbjct_taxonomy not in all_delta_taxo:
                    taxo_tree = insert_taxo_tot_delta(taxo_tree, query, pos_line, sbjct_taxonomy, args.rank, args.identical)
                    all_delta_taxo.append(sbjct_taxonomy)
        # print >>sys.stderr, p.get_memory_info()
        query_infos.pop(query)
    # print >>sys.stderr,  os.popen('ps -p %d -l' %(os.getpid())).read()
    # print >>sys.stderr,  int(os.popen('ps -p %d -o %s | tail -1' %(os.getpid(), "rss")).read())
    # print >>sys.stderr,  int(os.popen('ps -p %d -o %s | tail -1' %(os.getpid(), "vsz")).read())
    # print >>sys.stderr, p.get_memory_info()
    # args.tabfh.close()
    # ### output
    if args.tesxtfh:
        # print >>sys.stderr, 'beginning tree file writing', time.strftime("%y/%m/%d %H:%M:%S" , time.localtime(time.time()))
        try:
            tree_repr = rankoptimizerlib.to_tree(taxo_tree, query_name=False)
            print >>args.tesxtfh, tree_repr
        except IOError, err:
            print >>sys.stderr, err

    if args.jsonfh:
        try:
            kronaJson = rankoptimizerlib.KronaJSON(args.jsonfh, args.tabfh.name, taxo_tree)
            kronaJson.krona()
        except IOError, err:
            print >>sys.stderr, err

    if args.kronafh:
        # print >>sys.stderr, 'beginning xml krona file writing', time.strftime("%y/%m/%d %H:%M:%S" , time.localtime(time.time()))
        try:
            krona_xml = rankoptimizerlib.Krona(args.kronafh, args.tabfh.name, taxo_tree, args.kronahome, krona_js_on_server=args.jslocal)
            krona_xml.krona()
        except IOError, err:
            print >>sys.stderr, err

    if args.htmlxmlfh:
        # print >>sys.stderr, 'beginning html krona file writing', time.strftime("%y/%m/%d %H:%M:%S" , time.localtime(time.time()))
        try:
            krona_xml = rankoptimizerlib.Krona(args.htmlxmlfh, args.tabfh.name, taxo_tree, krona_url=args.kronahome, krona_js_on_server=args.jslocal)
            krona_xml.krona_html()
        except IOError, err:
            print >>sys.stderr, err

    if args.htmljsonfh:
        # print >>sys.stderr, 'beginning html krona file writing', time.strftime("%y/%m/%d %H:%M:%S" , time.localtime(time.time()))
        try:
            krona_json = rankoptimizerlib.KronaJSON(args.htmljsonfh, args.tabfh.name, taxo_tree, krona_url=args.kronahome, krona_js_on_server=args.jslocal)
            krona_json.krona_html()
        except IOError, err:
            print >>sys.stderr, err

    if args.dumpfh:
        import pickle
        pickle.dump(taxo_tree, args.dumpfh)
        args.dumpfh.close()

    if args.lca:
        print 'Lowest common ancestor:', lca(taxo_tree)
