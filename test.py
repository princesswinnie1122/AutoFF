import sys
import stk
import json
import re
import os
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QPushButton, QTableWidget, QTableWidgetItem, 
                             QLabel, QLineEdit, QTextEdit, QDialog, QFormLayout, QMessageBox, QHBoxLayout)
from PyQt6.QtCore import QSize
from cut_search import fragment_molecule, search_database
from link import link
from clean import clean_para, clean_all
from convert import info_to_dic
import py3Dmol
from rdkit import Chem


smiles = 'C1=CC=C2C(=C1)C3=CC=CC=C3N2C4=CC=C(C=C4)N(C5=CC=C(C=C5)N6C7=CC=CC=C7C8=CC=CC=C86)C9=CC=C(C=C9)N1C2=CC=CC=C2C2=CC=CC=C21'
self = fragment_molecule(smiles)
print(self)