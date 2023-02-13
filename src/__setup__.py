from os import path
from setuptools import setup, find_packages


this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), "r", encoding="utf-8") as f:
    long_description = f.read()

DESCRIPTION = (
    "Tools to annotate short reads taxonomy and function via pylogenetic tree evolutionary placement",
)
LONG_DESCRIPTION = (long_description,)
LONG_DESCRIPTION_CONTENT_TYPE = "text/markdown"
NAME = "phyloplacement"
AUTHOR = "Semidán Robaina Estévez, Jozé Manuel González Hernández"  # And other authors to be added
AUTHOR_EMAIL = "srobaina@ull.edu.es"
MAINTAINER = "Semidán Robaina Estévez"
MAINTAINER_EMAIL = "srobaina@gmail.com"
DOWNLOAD_URL = "http://github.com/robaina/TRAITS"
LICENSE = "Creative Commons Attribution 4.0 International"
VERSION = "0.0.1"

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    url=DOWNLOAD_URL,
    download_url=DOWNLOAD_URL,
    license=LICENSE,
    packages=find_packages(),
    install_requires=["numpy", "pyfastx"],  # Can't we do this wit envyronment.yml?
    entry_points={"console_scripts": ["phylolacement = phyloplacement.cli:main"]},
)
