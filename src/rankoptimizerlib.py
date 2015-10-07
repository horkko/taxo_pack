#! /usr/local/bin/python

# Corinne Maufrais
# Institut Pasteur, DSI/CIB
# maufrais@pasteur.fr
#


import os
import sys
try:
    LIB = os.environ['RANKOPTIMIZERLIB']
except:
    LIB = '/usr/local/bin'
    LIB = '../lib'

if LIB not in sys.path:
    sys.path.append(LIB)

try:
    SHARE = os.environ['RANKOPTIMIZERSHARE']
except:
    SHARE = '/Users/maufrais/Developpements2/taxo_pack/lib/KronaTools-2.1/src'

if SHARE not in sys.path:
    sys.path.append(SHARE)

global krona_js
krona_js = SHARE + '/krona-2.0.gensoft.js'

class Taxon(object):
    def __init__(self, name=None, rank=''):
        """
        t = Taxon(name)

        Creates a new instance of a binary tree with 'name object' (string).
        """
        self.name = name
        self.parent = None  # Taxon object
        self.childs = []    # list of Taxon object
        self.queriesS = []   # list of tuple [(query, pos_line)]
        self.repr = ''
        self.rank = rank
        self.nb_querys = 0

    def has_rank(self):
        return self.rank

    def has_child(self, child):
        """
        t.has_child(val) --> bool
        """
        for chld in self.childs:
            if chld.name == child:
                return True
        return False

    def has_childs(self):
        """
        t.has_childs() --> bool

        Return true if t has a list of childs, false otherwise. A child, it's a
        Taxon object.
        """
        return not self.childs == []

    def has_queries(self):
        """
        t.hasSpecies() --> bool
        """
        return self.queriesS != []

#     def has_query(self, name):
#         for qry, pos in self.queriesS:
#             if qry == name:
#                 return True
#         return False

#     def nb_query(self, name):
#         n = 0
#         for qry, pos in self.queriesS:
#             if qry == name:
#                 n += 1
#         return n

    def add_queries(self, query_all):
        self.queriesS.append(query_all)  # query_all = (query,  pos_line)

    def add_child(self, child, rank=''):
        """
        t.add_child(child)

        Creates a new Taxon object with child as name and append it as
        the childs list of t.
        """
        t = Taxon(child, rank)
        t.parent = self
        if self.childs == []:
            # t.repr = '|\t'
            t.repr = '|'
            self.childs = [t]
        else:
            # self.childs[-1].repr = '|\t'
            self.childs[-1].repr = '|'
            # t.repr = '.\t'
            t.repr = '.'
            self.childs.append(t)

    def get_child(self, child):
        for chld in self.childs:
            if chld.name == child:
                return chld

#     def nb_of_queries_in_childs(self):
#         nb = 0
#         for chld in self.childs:
#             nb += chld.nb_querys
#         return nb

    def has_one_child(self):
        return len(self.childs) == 1

# #################### Krona


class ElementXML(object):

    def __init__(self, outfh=None, indent=None):
        self.outfh = outfh
        self.indent = indent  # xml file indentation

    def start_elem(self, element, attr=''):
        # <element attributes>
        print >>self.outfh, '%s<%s%s>' % (self.space(), element, attr)

    def end_elem(self, element):
        # </element>
        print >>self.outfh, '%s</%s>' % (self.space(), element)

    def complete_elem(self, element, value, attr=''):
        # <element attributes>value</element>
        if value:
            print >>self.outfh, '%s<%s%s>%s</%s>' % (self.space(), element, attr, self.html_str(value), element)

#     def newLine(self):
#         print >>self.outfh, ''

    def increase_indent(self):
        self.indent += 1

    def decrease_indent(self):
        self.indent -= 1

    def get_indent(self):
        return self.indent

    def space(self):
        return '  ' * self.indent

    def html_str(self, value):
        if isinstance(value, str):
            # value = value.replace('&', '&amp;')
            value = value.replace('>', '&gt;')
            value = value.replace('<', '&lt;')
            return value
        else:
            return str(value)


class KronaDTD (ElementXML):
    def __init__(self, outfh=None, indent=None, krona_url='.', krona_js_on_server=False, collapse='true', key='true'):
        ElementXML.__init__(self, outfh=outfh, indent=indent)
        self.elems_attributes = {'krona': {'collapse': collapse, 'key': key},
                                 'attributes': {'magnitude': 'reads'},
                                 }
        self.krona_url = krona_url
        self.krona_js_on_server = krona_js_on_server

    def header_html_bad(self):
        print >>self.outfh, """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <meta charset="utf-8"/>
  <link rel="shortcut icon" href="%s/img/favicon.ico"/>
    <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"%s\"}</script>""" % (self.krona_url, self.krona_url)

        if not self.krona_js_on_server:
            print >>self.outfh, """<script src="%s/krona-2.0.js"></script>""" % self.krona_url
        else:
            print >>self.outfh, '<script type="text/javascript">'
            krona_jsfh = open(krona_js)
            line = krona_jsfh.readline()
            while line:
                print >>self.outfh, '    ', line.replace('\n', '')
                line = krona_jsfh.readline()
            print >>self.outfh, '</script>'
        print >>self.outfh, """</head>
<body>"""
        # if not self.krona_js_on_server:
        print >>self.outfh, """
  <img id="hiddenImage" src="%s/img/hidden.png" style="display:none"/>
  <img id="loadingImage" src="%s/img/loading.gif" style="display:none"/>""" % (self.krona_url, self.krona_url)
        print >>self.outfh, """
  <noscript>Javascript must be enabled to view this page.</noscript>
  <div style="display:none">"""

    def header_html(self, krona_jsfh=None):
        print >>self.outfh, """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <meta charset="utf-8"/>
  <link rel="shortcut icon" href="%s/img/favicon.ico"/>
    <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"%s\"}</script>""" % (self.krona_url, self.krona_url)

        if not self.krona_js_on_server:
            print >>self.outfh, """<script src="%s/krona-2.0.js"></script>""" % self.krona_url
        else:
            print >>self.outfh, '<script type="text/javascript">'
            line = krona_jsfh.readline()
            while line:
                print >>self.outfh, '    ', line.replace('\n', '')
                line = krona_jsfh.readline()
            print >>self.outfh, '</script>'
        print >>self.outfh, """</head>
<body>"""
        # if not self.krona_js_on_server:
        print >>self.outfh, """
  <img id="hiddenImage" src="%s/img/hidden.png" style="display:none"/>
  <img id="loadingImage" src="%s/img/loading.gif" style="display:none"/>""" % (self.krona_url, self.krona_url)
        print >>self.outfh, """
  <noscript>Javascript must be enabled to view this page.</noscript>
  <div style="display:none">"""


    def footer_html(self):
        print >>self.outfh, """
        </div></body></html>
        """

    def start_krona(self):
        # <krona ....>
        # krona={'collapse':'true', 'key':'true'}
        self.start_elem('krona', self._attr2str(self.elems_attributes['krona']))
        self.increase_indent()

    def end_krona(self):
        # </krona>
        self.decrease_indent()
        self.end_elem('krona')

    def color(self, color_values):
        # <color attribute="..." valueStart="..." valueEnd="..." hueStart="..." hueEnd="..." ...></color>
        # color_values = {attribute:'',  valueStart:'', valueEnd:'', hueStart:'',  hueEnd:''}
        self.complete_elem('color', '', self._attr2str(color_values))

    def datasets(self, datasets_values):
        # datasets_values:[set1, set2, ...]
        self.start_elem('datasets')
        self.increase_indent()
        for lds in datasets_values:
            # <dataset>set1</dataset>
            self.complete_elem('dataset', lds)
        self.decrease_indent()
        self.end_elem('datasets')

    def attributes(self, attributes_values):
        # 'attributes':{'magnitude': 'reads'},
        # <attributes magnitude="sample_attr_1">...</attributes>
        # attributes_values = {'sample_attr': [(sample_attr_name, attribute_values)],}
        #                      'sample_list': [(sample_list_name)],
        self.start_elem('attributes', self._attr2str(self.elems_attributes['attributes']))
        self.increase_indent()
        for san, av in attributes_values['sample_attr']:
            # <attribute attr >sample_attr_1</attribute>
            # <attribute display='' mono='' listNode='' listAll='' hrefBase='' target=''>...<attribute>
            # attribute_values = {display:'' mono:'' listNode:'' listAll:'' hrefBase:'' target:''}
            self.complete_elem('attribute', san, self._attr2str(av))
        for sln in attributes_values['sample_list']:
            # <list>sample_list_1</list>
            self.complete_elem('list', sln)
        self.decrease_indent()
        self.end_elem('attributes')

    def node(self, child_nodes):  # taxon
        # <node name='' href=''></node>
        # node_attributes={'name:'','href':''}
        node_attributes, rankNbRead_values, query_list_values, child_nodes = self._to_krona_node(child_nodes)
        self.start_elem('node', self._attr2str(node_attributes))
        self.increase_indent()
        # <sample_attr_1>
        # rankNbRead_values = [(name, [(v1,href),(v,href)...])] plusieurs val, si plusieurs dataset
        for san, sav in rankNbRead_values:
            self.sample_attr(san, sav)
        # <sample_list_1>

        for c in child_nodes:
            self.node(c)
        for sln, slv in query_list_values:
            self.sample_list(sln, slv)
        self.decrease_indent()
        self.end_elem('node')

    def _to_krona_node(self, taxon):
        node_attributes = {'name': taxon.name}
        if taxon.rank:
            rank_values = ('rank', [(taxon.rank, {})])
        else:
            rank_values = ('rank', [])

        child_nodes = taxon.childs
        # for c in taxon.childs:
        #    self.toKronaNode(c )

        read_values = ('reads', [(taxon.nb_querys, {})])
        # blast_values = ('blast',[(taxon.nb_querys, {})])
        list_values_qry = []
        # list_values_blast = []
        if taxon.has_queries():
            for query in taxon.queriesS:
                list_values_qry.append(("%-40s\t%s" % (query[0], query[1]), {}))
                # list_values_blast.append((query[1], {}))
        # list_values = [('read_members',list_values_qry), ('blast_members',list_values_blast)]
        query_list_values = [('read_members', list_values_qry)]

        # return node_attributes, [read_values, blast_values, rank_values], list_values, child_nodes
        return node_attributes, [read_values, rank_values], query_list_values, child_nodes

    def sample_attr(self, sample_attr_name, sample_attr_values):
        if sample_attr_values:
            self.start_elem(sample_attr_name)
            self.increase_indent()
            for vl, av in sample_attr_values:  # plusieurs val, si plusieurs dataset
                # <val href=''>7</val>
                self.val(vl, av)
            self.decrease_indent()
            self.end_elem(sample_attr_name)

    def val(self, value, val_attribute):
        # <val href=''>12</val>
        self.complete_elem('val', value, self._attr2str(val_attribute))

    def sample_list(self, sple_list_name, sple_list_values):
        if sple_list_values:
            self.start_elem(sple_list_name)
            self.increase_indent()
            self.start_elem('vals')
            self.increase_indent()
            for vl, av in sple_list_values:
                # <val herf=''>7</val>
                self.val(vl, av)
            self.decrease_indent()
            self.end_elem('vals')
            self.decrease_indent()
            self.end_elem(sple_list_name)

    def _attr2str(self, attrs):
        s = ''
        for k, v in attrs.items():
            if v:
                s += ' %s="%s"' % (k, str(v))
        return s

    def _setAttr(self, elem, attr, value):
        self.attributes[elem][attr] = value


class Krona (KronaDTD):
    def __init__(self,  outfh=None, file_name=None, taxo_tree=None, krona_url='.', krona_js_on_server=False, collapse='false', key='true'):
        """
        * Object for translating one treeobject into one Krona xml file.
        """
        KronaDTD.__init__(self, outfh=outfh, indent=0, krona_url=krona_url, krona_js_on_server=krona_js_on_server, collapse=collapse, key=key)
        self.taxo_tree = taxo_tree
        self.file_name = file_name

    def krona(self):
        self.indent = 0
        # xmlKrona = KronaDTD(self.outfh, self.indent)
        self.start_krona()
        self.datasets([self.file_name])
        # color_values = {'attribute':'', 'default':'true'}
        # xmlKrona.color(color_values)
        attributes_values = {'sample_attr': [('reads', {'display': 'Nb of reads', 'listAll': 'read_members'}),
                                             # ('blast', {'display':'Blast offset', 'listAll':'blast_members'}),
                                             ('rank', {'display': 'Rank', 'mono': 'true'})],
                             # 'sample_list':[]
                             'sample_list': ['read_members',
                                             'blast_members']
                             }

        self.attributes(attributes_values)

        self.node(self.taxo_tree)
        self.end_krona()

    def krona_html(self, krona_jsfh=None):
        self.indent = 0
        # xmlKrona = KronaDTD(self.outfh, self.indent, self.krona_url, self.krona_js_on_server)
        self.header_html(krona_jsfh)
        self.start_krona()
        self.datasets([self.file_name])
        # color_values = {'attribute':'', 'default':'true'}
        # xmlKrona.color(color_values)
        attributes_values = {'sample_attr': [('reads', {'display': 'Nb of reads', 'listAll': 'read_members'}),
                                             # ('blast', {'display':'Blast offset', 'listAll':'blast_members'}),
                                             ('rank', {'display': 'Rank', 'mono': 'true'})],
                             # 'sample_list':[]
                             'sample_list': ['read_members',
                                             'blast_members']
                             }

        self.attributes(attributes_values)
        self.node(self.taxo_tree)
        self.end_krona()
        self.footer_html()
    # ####################


class ElementJSON(object):
    """
    XML                              JSON
    <e/>                             "e": null
    <e>text</e>                      "e": "text"
    <e attr="value" />               "e": {"@attr": "value"}    @ ou _
    <e attr="value">text</e>         "e": {"@attr": "value", "#text": "text"}
    <e> <a>text</a ><b>text</b> </e> "e": {"a": "text", "b": "text"}
    <e> <a>text</a> <a>text</a> </e> "e": {"a": ["text", "text"]}
    <e> text <a>text</a> </e>        "e": {"#text": "text", "a": "text"}
    """

    def __init__(self, outfh=None, indent=None):
        self.outfh = outfh
        self.indent = indent  # xml file indentation

    def space(self):
        return ''

    def start_brace(self):
        # {
        print >>self.outfh, '{',

    def end_brace(self, coma=''):
        # }
        print >>self.outfh, '}%s' % (coma),

    def start_element(self, element):
        # element:
        print >>self.outfh, '"%s":' % (element),

#     def end_element(self, coma=''):
#         # element:
#         print >>self.outfh, '%s' % (coma),

    def elem_attributes(self, k_attr={}, coma=''):
        # krona = {'collapse':'true', 'key':'true'}
        # 'collapse':'true', 'key':'true',
        st = ''
        for k, v in k_attr.items():
            st += '"%s":"%s",' % (k, v)
        print >>self.outfh, st[:-1] + coma,

    def complete_simple_element(self, element, value, coma=''):
        # <e>text</e>  ==>   "e": "text"
        print >>self.outfh, '"%s":%s%s' % (element, value, coma),

    def complete_simple_element_str(self, element, value, coma=''):
        # <e>text</e>  ==>   "e": "text"
        return '"%s":%s%s' % (element, value, coma)

    def complete_simple_dbl_quote_element(self, element, value, coma=''):
        # <e>text</e>  ==>   "e": "text"
        print >>self.outfh, '"%s":"%s"%s' % (element, value, coma),

    def complete_simple_dbl_quote_element_str(self, element, value, coma=''):
        # <e>text</e>  ==>   "e": "text"
        return '"%s":"%s"%s' % (element, value, coma)

    def complete_complex_brace_element(self, element, value_str, coma=''):
        # <e>text</e>  ==>   "e": "text"
        print >>self.outfh, '"%s":{%s}%s' % (element, value_str, coma),

    def complete_complex_brace_element_str(self, element, value_str, coma=''):
        # <e>text</e>  ==>   "e": "text"
        return '"%s":{%s}%s' % (element, value_str, coma),

    def complete_element_with_attr(self, element, value, attrs={}, coma=''):
        # <e attr="value">text</e>   ==>   "e": {"attr": "value", "#text": "text"}
        print >>self.outfh, '"%s":{' % (element),
        self.elem_attributes(self, attrs)
        print >>self.outfh, '"__text":"%s"' % (value),
        print >>self.outfh, '}%s' % (coma),

    def start_list_element(self, element):
        # element: [
        print >>self.outfh, '"%s":[' % (element),

    def end_list_element(self, coma=''):
        print >>self.outfh, ']%s' % (coma),

    def increase_indent(self):
        self.indent += 1

    def decrease_indent(self):
        self.indent -= 1

    def get_indent(self):
        return self.indent


class KronaJSONDTD (ElementJSON):

    def __init__(self, outfh=None, indent=None, krona_url='.', krona_js_on_server=False, collapse='true', key='true'):
        ElementJSON.__init__(self, outfh=outfh, indent=indent)
        self.elems_attributes = {'krona': {'_collapse': collapse, '_key': key},
                                 'attributes': {'_magnitude': 'reads'},
                                 }

        self.krona_url = krona_url
        self.krona_js_on_server = krona_js_on_server

    def header_html_bad(self):
        print >>self.outfh, """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <meta charset="utf-8"/>
  <link rel="shortcut icon" href="%s/img/favicon.ico"/>
    <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"%s\"}</script>""" % (self.krona_url, self.krona_url)
        print >>self.outfh, """<script type='text/javascript' src='https://code.jquery.com/jquery-git.js'></script>
    <script type='text/javascript' src='https://x2js.googlecode.com/hg/xml2json.js'></script>"""
        if not self.krona_js_on_server:
            print >>self.outfh, """<script src="%s/krona-2.0.js"></script>""" % self.krona_url
        else:
            print >>self.outfh, '<script type="text/javascript">'
            krona_jsfh = open(krona_js)
            line = krona_jsfh.readline()
            while line:
                print >>self.outfh, '    ', line.replace('\n', '')
                line = krona_jsfh.readline()
        print >>self.outfh, '<script>'
        print >>self.outfh, "var jsontest = "


    def header_html(self, krona_jsfh=None):
        print >>self.outfh, """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <meta charset="utf-8"/>
  <link rel="shortcut icon" href="%s/img/favicon.ico"/>
    <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"%s\"}</script>""" % (self.krona_url, self.krona_url)
        print >>self.outfh, """<script type='text/javascript' src='https://code.jquery.com/jquery-git.js'></script>
    <script type='text/javascript' src='https://x2js.googlecode.com/hg/xml2json.js'></script>"""
        if not self.krona_js_on_server:
            print >>self.outfh, """<script src="%s/krona-2.0.js"></script>""" % self.krona_url
        else:
            print >>self.outfh, '<script type="text/javascript">'
            line = krona_jsfh.readline()
            while line:
                print >>self.outfh, '    ', line.replace('\n', '')
                line = krona_jsfh.readline()
        print >>self.outfh, '<script>'
        print >>self.outfh, "var jsontest = "


    def header_html2(self):
        print >>self.outfh, ';'
        print >>self.outfh, """
        var x2js = new X2JS();
    jsontest = JSON.stringify(jsontest);
    var j = JSON.parse(jsontest);
    <!-- chargement dans le DOM pour l'utiliseration par le js de krona //-->
    var xmlDoc = x2js.json2xml_str(j);
        """
        print >>self.outfh, '</script>'
        print >>self.outfh, """</head>
<body>"""

        # if not self.krona_js_on_server:
        print >>self.outfh, """
  <img id="hiddenImage" src="%s/hidden.png" style="display:none"/>
  <img id="loadingImage" src="%s/loading.gif" style="display:none"/>""" % (self.krona_url, self.krona_url)
        print >>self.outfh, """
  <noscript>Javascript must be enabled to view this page.</noscript>
  <div style="display:none">
  <div id='krona' style="display:none">
    <?xml version="1.0" encoding="UTF-8" ?>
    </div>

    <script>
        var kronacontent = $("#krona");
        kronacontent.append(xmlDoc);
    </script>
  """

    def footer_html(self):
        print >>self.outfh, """
        </div></body></html>
        """

    def start_krona_with_attr(self):
        # <krona collapse='true' key = 'true'...> ==> {krona : {'collapse':'true', 'key':'true',
        # krona={'collapse':'true', 'key':'true'}
        self.start_brace()
        self.start_element('krona')
        self.start_brace()
        if 'krona' in self.elems_attributes:
            self.elem_attributes(self.elems_attributes['krona'], coma=',')

    def end_krona(self):
        # </krona> ==>}}
        self.end_brace()
        self.end_brace()

    def datasets(self, datasets_values):
        # datasets_values:[set1, set2, ...]
        # <datasets><dataset>text1</dataset><dataset>text2</datatset></datasets>
        # <e> <a>text1</a> <a>text2</a> </e> ==> "e": {"a": ["text1", "text2"]}
        self.complete_simple_dbl_quote_element('dataset', datasets_values, coma=',')

    def attributes(self):
        # <attributes magnitude="reads"> ... </attributes> ==> attributes: {magnitude: "reads" ,
        # 'attributes':{'magnitude': 'reads'},
        self.start_element('attributes')
        self.start_brace()
        self.elem_attributes(self.elems_attributes['attributes'], coma=',')

        # <attribute display="Nb of reads" listAll="read_members">reads</attribute> ==> attribute: [{display:"Nb of reads" , listAll:"read_members", #text:"reads"},
        # <attribute mono="true" display="Rank">rank</attribute>                    ==>   {"mono":"true", "display":"Rank", #text : "rank"}],
        self.start_list_element('attribute')
        # ('reads', {'display':'Nb of reads', 'listAll':'read_members'})
        self.start_brace()
        self.elem_attributes({'_display': 'Nb of reads', '_listAll': 'read_members'}, coma=',')
        print >>self.outfh, '"__text":"%s"' % ('reads'),
        self.end_brace(coma=',')

        # ('rank', {'display':'Rank', 'mono':'true'})
        self.start_brace()
        self.elem_attributes({'_display': 'Rank', '_mono': 'true'}, coma=',')
        print >>self.outfh, '"__text":"%s"' % ('rank'),
        self.end_brace()

        self.end_list_element(coma=',')

        # <list>read_members</list> ==> list:  read_members
        self.complete_simple_dbl_quote_element('list', 'read_members')

        # </attributes> ==>}}
        self.end_brace(coma=',')

    def node(self, node, coma=','):  # taxon
        node_attributes, nb_reads, rank, read_members, child_nodes = self._to_krona_node(node)
        # node_attributes={'name:'','href':''}
        # <node name='' href=''> ... </node> ==> ==> node: {name: "" , href: "",

        self.start_list_element('node')

        self.start_brace()
        self.elem_attributes(node_attributes, ',')

        # <sample_attr_1> ==> sample_attr_1:

        self.complete_complex_brace_element('reads',  self.complete_simple_dbl_quote_element_str('val',  nb_reads, ''), ',')
        if rank:
            self.complete_complex_brace_element('rank',  self.complete_simple_dbl_quote_element_str('val',  rank, ''), ',')

        if len(child_nodes) > 1:
            self.start_list_element('node')
            for c in child_nodes[:-1]:
                self.node_without_node(c, ',')
                print >>self.outfh, ',',
            self.node_without_node(child_nodes[-1], '')
            self.end_list_element()
        elif child_nodes:
            if read_members:
                coma = ','
            else:
                coma = ','  # coma =''
            self.node(child_nodes[0], coma)
            # coma = ''
        # <sample_list_1>
        # read_members = {'vals':{'vals':[a,b,c,d]}}

#        if read_members:
#            self.start_list_element('read_members')
#            for info_dict in read_members[:-1]:
#                self.elem_attributes(info_dict, '')
#            self.elem_attributes(read_members[0])
#            self.end_list_element()

        if read_members:
            self.start_element('read_members')
            self.start_brace()
            self.complete_complex_brace_element('vals',  self.complete_simple_element_str('val',  read_members, ''), '')
            self.end_brace(coma=',')

        self.end_brace()
        self.end_list_element(coma)

    def node_without_node(self, node,  coma=','):  # taxon
        node_attributes, nb_reads, rank, read_members, child_nodes = self._to_krona_node(node)
        # node_attributes={'name:'','href':''}
        # <node name='' href=''> ... </node> ==> ==> node: {name: "" , href: "",

        self.start_brace()
        self.elem_attributes(node_attributes, ',')

        # <sample_attr_1> ==> sample_attr_1:
        self.complete_complex_brace_element('reads',  self.complete_simple_dbl_quote_element_str('val',  nb_reads, ''), ',')
        if rank:
            self.complete_complex_brace_element('rank',  self.complete_simple_dbl_quote_element_str('val',  rank, ''), ',')

        if len(child_nodes) > 1:
            self.start_list_element('node')
            for c in child_nodes[:-1]:
                self.node_without_node(c, ',')
                print >>self.outfh, ',',
            self.node_without_node(child_nodes[-1], '')
            self.end_list_element()
        elif child_nodes:
            self.node(child_nodes[0], '')
        # <sample_list_1>
        # read_members = {'vals':{'vals':[a,b,c,d]}}
        if read_members:
            self.start_element('read_members')
            self.start_brace()
            self.complete_complex_brace_element('vals',  self.complete_simple_element_str('val',  read_members, ''), '')
            self.end_brace(coma=',')

        self.end_brace()

    def _to_krona_node(self, taxon):
        node_attributes = {'_name': taxon.name.replace(':', '_').replace('"', '_')}
        if taxon.rank:
            rank_values = taxon.rank
        else:
            rank_values = ''

        child_nodes = taxon.childs
        # for c in taxon.childs:
        #    self.toKronaNode(c )

        list_values_qry = []
        if taxon.has_queries():
            for query in taxon.queriesS:
                list_values_qry.append(query[0].replace(':', '_').replace('"', '_') + '\t%s' % query[1])
                # list_values_qry.append({'name':query[0].replace(':', '_').replace('"', '_') + '\t%s' % query[1]})
        read_members = list_values_qry

        return node_attributes, taxon.nb_querys, rank_values, read_members, child_nodes


class KronaJSON(KronaJSONDTD):
    def __init__(self,  outfh=None, file_name=None, taxo_tree=None, krona_url='.', krona_js_on_server=False,  collapse='false', key='true'):
        """
        * Object for translating one treeobject into one Krona xml file.
        """
        KronaJSONDTD.__init__(self, outfh=outfh, indent=0, krona_url=krona_url, krona_js_on_server=krona_js_on_server, collapse=collapse, key=key)
        self.taxo_tree = taxo_tree
        self.file_name = file_name

    def krona(self):
        self.indent = 0
        self.start_krona_with_attr()
        self.datasets(self.file_name)
        self.attributes()
        self.node(self.taxo_tree)
        self.end_krona()

    def krona_html(self, krona_jsfh=None):
        self.indent = 0
        self.header_html(krona_jsfh)
        self.indent = 0
        self.start_krona_with_attr()
        self.datasets(self.file_name)
        self.attributes()
        self.node(self.taxo_tree)
        self.end_krona()
        self.header_html2()
        self.footer_html()

#import newickTree_rankopti 


# def toNewickTree(taxon, tree, Acc):
#     p = 0
#     # print len(taxon.childs)
#     for taxa in taxon.childs:
#         if not tree.left:
#             tree.insertLeft(None)
#             # tree.left.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '') + '#' + taxa.rank + ':' + str(len(taxa.queriesS))
#             tree.left.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '') + '#' + taxa.rank + ':' + str(taxa.nb_querys)
#             tree.left = toNewickTree(taxa, tree.left, Acc)
#         elif not tree.right:
#             tree.insertRight(None)
#             # tree.right.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '') + '#' + taxa.rank + ':' + str(len(taxa.queriesS))
#             tree.right.root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '') + '#' + taxa.rank + ':' + str(taxa.nb_querys)
#             tree.right = toNewickTree(taxa, tree.right, Acc)
#         else:
#             p = tree.appendOther(None)
#             # tree.other[p].root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '') + '#' + taxa.rank + ':' + str(len(taxa.queriesS))
#             tree.other[p].root = taxa.name.replace(' ', '_').replace(':', '_').replace('.', '') + '#' + taxa.rank + ':' + str(taxa.nb_querys)
#             tree.other[p] = toNewickTree(taxa, tree.other[p], Acc)
#             p += 1
#     return tree
# 
# 
# def toDnd(taxon, Acc=False):
#     tree = newickTree_rankopti.NewickTree('root')
#     return toNewickTree(taxon, tree, Acc)


def _to_tree(taxon, name=True, s='', sv='', sd='.', query_name=False):
    if name:
        if taxon.rank:
            s += '+ ' + taxon.name + ' (' + taxon.rank + ')' + ';'
        else:
            s += '+ ' + taxon.name + ';'

    if len(taxon.childs) > 1:
        s += '#' + str(taxon.nb_querys) + '\n'
        for c in taxon.childs:
            s, sv = _to_tree(c, True, s + sd, sv, sd + c.repr, query_name)

    elif len(taxon.childs) == 1:
        c = taxon.childs[0]
        if c.rank:
            s += c.name + ' (' + c.rank + ')' + ';'
        else:
            s += c.name + ';'
        s, sv = _to_tree(c, False, s, sv, sd, query_name)
    else:
        # s += '#' + str(len(taxon.queriesS)) + '\n'
        s += '#' + str(taxon.nb_querys) + '\n'
        if taxon.has_queries():
            if query_name:
                for query, pos_line in taxon.queriesS:
                    s += "%s\n" % ((sd + ' - ' + query + ' (' + str(pos_line) + ') '))
            # else:
            #    s += "%s\n" % ((sd + ' - ' +  str(len(taxon.queriesS)) + ' queries'))

    return s, sv


def to_tree(taxon, name=True, s='', sv='', sd='.', query_name=False):
    s, sv = _to_tree(taxon, name, s, sv, sd, query_name)
    s += sv + '\n'
    return s
