# Author: Emmanule Quevillon, tuco@pasteur.fr
# Institut Pasteur, Centre d'informatique pour les biologistes
from __future__ import print_function
import sys


class Error(object):
    def __init__(self, err):
        self.err = str(err)

    def __repr__(self):
        print("[ERROR] %s " % self.err, file=sys.stderr)

    @staticmethod
    def fatal(msg):
        """
        Prints a FATAL message on STDERR and exit with code 1
        :param msg: Message to print
        :type msg: str
        :return:
        """
        print("[FATAL] %s" % str(msg), file=sys.stderr)
        sys.exit(1)


class GoldenError(object):
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[GoldenError] " + self.err


class TaxOptimizerError(object):
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[TaxOptimizerError] " + self.err


class ParserError(object):
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[ParserError] " + self.err


class ParserWarning:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[ParserWarning] " + self.err
