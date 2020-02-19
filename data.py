"""
Tool that contains the utils to export and import Pickle files
"""

import pickle
import logging
import os

__author__ = "BATIER Lucas, GONZALEZ JIMENEZ Alvaro, PILIPOVIC Predrag"
__copyright__ = "Copyright 2020, The Salmonela Project"
__license__ = "GPL"
__version__ = "1.0.1"
__email__ = "lucas.batier@hotmail.fr, alvaro.gonzalez-jimenez@grenoble-inp.org, predrag.pilipovic@grenoble-inp.org"
__status__ = "Production"


def save_obj(obj, name):
    """
    @param obj: Dictionary of k-mers
    @param name: Name of the output file
    @return: Create a pickle file with the content the obj
    """
    dirpath = os.getcwd()
    with open(dirpath + "/" + name, "wb") as f:
        logging.info("Saving pickle file with name: " + str(name) + ".pkl")
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    """
    @param name: Name of the pickle file to read
    @return: Load in a dictionary the content from the pickle file
    """
    dirpath = os.getcwd()
    with open(dirpath + "/" + name, "rb") as f:
        logging.info("Loading pickle file with name: " + str(name) + ".pkl")
        return pickle.load(f)
