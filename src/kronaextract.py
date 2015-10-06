#! /usr/local/bin/python

# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# maufrais@pasteur.fr
#

# version 2.0

import xml.etree.ElementTree as ET
import sys
import argparse


def extract_reads_from_all_children(nodes, list_of_reads):
    for node in nodes:
        # reads = node.find('reads')
        # nb_reads = int(reads.find('val').text)
        read_m = node.find('read_members')
        if read_m is not None:
            vals = read_m.find('vals')
            # all_vals = vals.findall('val')
            list_of_reads.extend(vals.findall('val'))
        newnodes = node.findall('node')
        list_of_reads = extract_reads_from_all_children(newnodes, list_of_reads)
    return list_of_reads


if __name__ == '__main__':

    usage = "kronaextract [options] -i <FILE>  -n <STRING>"
    epilog = """
    kronaextract extract list of reads and blast offset for a given taxonomic name, from a xml file obtained by the rankoptimizer program.
The output file could be split into two files: the 'prefix.seq' file contains reads names and the 'prefix.offset' file contains the corresponding taxoptimizer's line offset.

kronaextract use Krona 2.1, an interactive metagenomic visualization tool in a Web browser.  (http://sourceforge.net/p/krona/home/krona/):
    Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.
"""
    parser = argparse.ArgumentParser(prog='kronaextract.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage, epilog=epilog)
    general_options = parser.add_argument_group(title="Options", description=None)
    general_options.add_argument("-i", "--in",
                                 dest="krona_xml_file",
                                 help="Xml input file with Krona 2.1 Specification. A rankoptimizer ouptut file is recommended)",
                                 metavar="file",
                                 required=True)
    general_options.add_argument("-n", "--taxo_name",
                                 dest="taxoname", metavar="STRING",
                                 help="Taxonomic name.")
    general_options.add_argument("-o", "--out", dest="outfile",
                                 help="Output file.",
                                 type=argparse.FileType('w'),
                                 metavar="file")
    general_options.add_argument("-s", "--split_prefix",
                                 dest="prefix", metavar="str",
                                 help="Split output file into two files with the given prefix name")

    args = parser.parse_args()

    bl_line = 0

    # ===== Tabulated file parsing
    try:
        xmltree = ET.parse(args.krona_xml_file)
    except IOError, err:
        print >>sys.stderr, err
        sys.exit(0)

    root = xmltree.getroot()

    nodes = root.findall(".//node/[@name='%s']" % args.taxoname)
    list_of_reads = []
    if not nodes:
        print >>sys.stderr, 'No result for: %s' % args.taxoname
        sys.exit()
    for nd in nodes:
        reads = nd.find('reads')
        nb_reads = int(reads.find('val').text)
        read_m = nd.find('read_members')
        if read_m is not None:
            vals = read_m.find('vals')
            list_of_reads = vals.findall('val')

        nodes = nd.findall('node')
        list_of_reads = extract_reads_from_all_children(nodes, list_of_reads)
#        if len(list_of_reads) !=  nb_reads:
#                print 'Un probleme  a regler: %s reads reel, %s reads theorique' % (len(list_of_reads), nb_reads)
#        else:
#            print 'tout est ok: %s reads' % len(list_of_reads)

    if args.prefix:
        outfh_name = open(args.prefix + '.seq', 'w')
        outfh_offset = open(args.prefix + '.offset', 'w')

    for read_info in list_of_reads:
        if args.outfile:
            print >>args.outfile, read_info.text
        if args.prefix:
            fld = read_info.text.split('\t')
            print >>outfh_name, fld[0]
            print >>outfh_offset, fld[1]
