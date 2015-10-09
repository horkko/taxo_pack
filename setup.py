#!/usr/bin/env python
# -*- coding: utf-8 -*-
from distutils.core import setup
from distutils.command.install_egg_info import install_egg_info
from distutils.versionpredicate import VersionPredicate
from distutils.command.build import build
import sys


class nohup_egg_info(install_egg_info):
    def run(self):
        # there is nothing to install in sites-package
        # so I don't put any eggs in it
        pass


class check_and_build(build):
    def run(self):
        chk = True
        for req in require_pyt:
            chk &= self.chkpython(req)
        for req in require_mod:
            chk &= self.chkmodule(req)
        if not chk:
            sys.exit(1)
        build.run(self)

    def chkpython(self, req):
        chk = VersionPredicate(req)
        ver = '.'.join([str(v) for v in sys.version_info[:2]])
        if not chk.satisfied_by(ver):
            print >> sys.stderr, "Invalid python version, expected %s" % req
            return False
        return True

    def chkmodule(self, req):
        chk = VersionPredicate(req)
        try:
            mod = __import__(chk.name)
        except:
            print >> sys.stderr, "Missing mandatory %s python module" % chk.name
            return False
        for v in ['__version__', 'version']:
            ver = getattr(mod, v, None)
            break
        try:
            if ver and not chk.satisfied_by(ver):
                print >> sys.stderr, "Invalid module version, expected %s" % req
                return False
        except:
            pass
        return True

require_pyt = ['python (>=2.7, <3.0)']
require_mod = []


setup(name="taxo_pack",
      version="2.0",
      author="Corinne Maufrais",
      author_email="corinne.maufrais@pasteur.fr",
      license='GPLv3',
      description=("""All programs permit to analyze and visualize the taxonomic abundance
of set of sequences compared against a sequence database with BLAST.
Taxoptimizer parse the blast output report and add the NCBI Taxonomy database information to each HSP.
Rankoptimizer analyze the taxonomy abundance of a set of sequences, pre-process by the taxoptimizer program, and format result with Krona,
an interactive metagenomic visualization in a Web browser.
kronaextract extract sub-list of Query ID from a set of sequences matching a given taxon name and/or their offset number in the blast report """),
      scripts=['src/taxoptimizer', 'src/taxoptimizer', 'src/kronaextract'],
      cmdclass={'install_egg_info': nohup_egg_info},
      package_dir={'': 'src'},
      install_requires=['golden >= 3.0',
                        'bsddb3 >=6.1.0',
                        'KronaTools >= 2.6',
                        ],
      )
