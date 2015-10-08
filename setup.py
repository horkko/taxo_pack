try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name = "taxo_pack",
    version = "2.0",
    description = (""),
    license='',
    author = "Corinne Maufrais",
    author_email = "corinne.maufrais@pasteur.fr",
    package_dir={'': 'src'},
    install_requires = ['golden >= 3.0', 
	'bsddb3 >=6.1.0', 
	'KronaTools >= 2.6', 	
],
) 
