"""
GUI app to interact with the models easily
"""
import sys
from PIL.ImageQt import ImageQt
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QLabel, QWidget, QVBoxLayout, QLineEdit, QTabWidget, QPushButton, QFileDialog
from PyQt5.QtGui import QImage, QPixmap
from rdkit import Chem
from math import floor, sqrt
from rdkit.Chem import Draw, rdMolDescriptors, Descriptors
from chemical_models import *
from graphs.molecule_relative_similarity import mol_similarity_grid

# Load models
logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')
melting_point_model = MeltingPoint('melting_gse')


class MainWidget(QTabWidget):
    def __init__(self):
        QTabWidget.__init__(self)

        # Stores the list of SMILES strings
        self.mols = []

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.setWindowTitle('Chemistry')

        self.selection_instruction = QLabel()
        self.select_smiles_file = QPushButton('Select File')
        self.mol_img = QLabel()
        self.logP = QLabel()
        self.solubility = QLabel()
        self.melting = QLabel()

        self.select_smiles_file.pressed.connect(self.select_file)

        self.smiles = QLineEdit()
        self.smiles.setPlaceholderText('Enter SMILES string separated by commas')
        self.smiles.returnPressed.connect(self.string_returned)

        self.mol_selection = QWidget()
        self.mol_display = QWidget()
        self.mol_graph = QWidget()

        self.addTab(self.mol_selection, 'Select')
        self.addTab(self.mol_display, 'View')
        self.addTab(self.mol_graph, 'Graph')

        self.mol_selection_ui()
        self.mol_display_ui()

    def mol_selection_ui(self):
        layout = QVBoxLayout()
        layout.addWidget(self.smiles)
        layout.addWidget(self.selection_instruction)
        layout.addWidget(self.select_smiles_file)
        self.selection_instruction.setText('Or select a file containing a list of '
                                           'SMILES strings separated by commas')
        self.mol_selection.setLayout(layout)

    def mol_display_ui(self):
        layout = QVBoxLayout()
        layout.addWidget(self.mol_img)
        self.mol_display.setLayout(layout)

    def mols_loaded(self):
        if len(self.mols) == 1:
            img = QImage(ImageQt(Draw.MolToImage(self.mols[0])))
            self.mol_img.setPixmap(QPixmap(img))
        elif len(self.mols) <= 35:
            fps = []
            for mol in self.mols:
                fps.append((rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol), mol))

            dim = int(int(len(self.mols)) / floor(sqrt(len(self.mols))))
            img = QImage(ImageQt(mol_similarity_grid(fps, (int(700/dim), int(700/dim)), dim)))
            self.mol_img.setPixmap(QPixmap(img))

    @pyqtSlot()
    def string_returned(self):
        self.mols = []
        for mol in self.smiles.text().replace(' ', '').split(','):
            try:
                self.mols.append(Chem.MolFromSmiles(mol))
            except Exception as e:
                print(e)
        self.mols_loaded()

    def select_file(self):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        filenames = []

        if dlg.exec_():
            filenames = dlg.selectedFiles()
            f = open(filenames[0], 'r')

            with f:
                data = f.read()
                self.contents.setText(data)

        # img = QImage(ImageQt(Draw.MolToImage(mol)))
        # fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
        # logP = logP_model.run(fp)
        # logP_sol = logP_solubility_model.run(logP)
        # atom_pair_sol = atom_pair_sol_model.run(fp)
        # combined_sol = combined_model.run(mol, logP, logP_sol, atom_pair_sol)
        # mg_l_sol = (10.0**float(combined_sol)) * float(Descriptors.ExactMolWt(mol))
        # mp = melting_point_model.run(combined_sol, logP)
        #
        # self.logP.setText('LogP: ' + str(round(logP, 2)))
        # self.solubility.setText('Water Solubility(mg/L): ' + str(round(mg_l_sol, 2)))
        # self.melting.setText('Melting Point(C): ' + str(round(mp, 2)))
        # self.mol_img.setPixmap(QPixmap(img))


app = QApplication(sys.argv)

widget = MainWidget()
widget.resize(800, 800)
widget.show()
sys.exit(app.exec_())
