# GUI app to interact with the models easily
import sys
from PIL.ImageQt import ImageQt
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QImage, QPixmap
from math import sqrt, ceil
from rdkit.Chem import Draw
from chemical_models import *
from functions import logs_to_mg_ml
from graphs.molecule_relative_similarity import mol_similarity_grid

logP_model = LogP('logP')
logP_solubility_model = LogPSolubility('logS_logP')
atom_pair_sol_model = AtomPairSolubility('water_solubility')
combined_model = CombinedSolubility('combined_solubility')
melting_point_model = MeltingPoint('melting_gse')

w, h = 800, 600


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Chemistry')
        self.resize(w, h)

        self.central_widget = QWidget()
        self.layout = QVBoxLayout()
        self.central_widget.setLayout(self.layout)

        # Stores both the RDKit Mol and SMILES
        self.mols = []
        self.smiles = []

        self.toolbar = QToolBar()

        self.mol_img = QLabel()
        self.mol_name = QLabel()
        self.logP = QLabel()
        self.solubility = QLabel()
        self.melting = QLabel()

        self.smiles_entry = EnterSmiles(self.mols, self.smiles, self.mols_loaded)

        self.create_toolbar()

        self.layout.addWidget(self.mol_img)
        self.layout.addWidget(self.mol_name)
        self.layout.addWidget(self.logP)
        self.layout.addWidget(self.solubility)
        self.layout.addWidget(self.melting)

        self.setCentralWidget(self.central_widget)

        self.mols_loaded()

    def create_toolbar(self):
        menu_bar = self.menuBar()
        menu_bar.setNativeMenuBar(True)

        file_menu = menu_bar.addMenu(' &File')

        smiles_act = QAction(' Enter SMILES', self)

        smiles_act.setShortcut('Ctrl+J')
        smiles_act.setStatusTip('Manually enter SMILES strings')
        smiles_act.triggered.connect(self.smiles_action)

        file_menu.addAction(smiles_act)
        file_menu.addAction('Open...')

        self.setMenuBar(menu_bar)

    def mols_loaded(self):
        if len(self.mols) == 0:
            self.mol_img.setText('No molecules selected')
        elif len(self.mols) == 1:
            img = QImage(ImageQt(Draw.MolToImage(self.mols[0], size=(int(w/2), int(h/2)))))
            self.mol_img.setPixmap(QPixmap(img))

            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(self.mols[0])
            logP = logP_model.run(fp)
            logP_sol = logP_solubility_model.run(logP)
            atom_pair_sol = atom_pair_sol_model.run(fp)
            combined_sol = combined_model.run(self.mols[0], logP,
                                              logP_sol, atom_pair_sol)
            mg_ml_sol = logs_to_mg_ml(combined_sol, self.mols[0])
            mp = melting_point_model.run(combined_sol, logP)

            self.mol_name.show()
            self.logP.show()
            self.solubility.show()
            self.melting.show()

            self.mol_name.setText(self.smiles[0])
            self.logP.setText('LogP: ' + str(round(logP, 2)))
            self.solubility.setText('Water Solubility(mg/mL): ' +
                                    str(round(mg_ml_sol, 2)))
            self.melting.setText('Melting Point(C): ' + str(round(mp, 2)))
            self.mol_img.setPixmap(QPixmap(img))

        elif len(self.mols) <= 35:
            fps = []
            for mol in self.mols:
                fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
                fps.append((fingerprint, mol))
            dim = int(ceil(sqrt(len(self.mols))))

            img = QImage(ImageQt(mol_similarity_grid(fps, (int(700/dim),
                                                     int(700/dim)), dim)))

            self.mol_name.hide()
            self.logP.hide()
            self.solubility.hide()
            self.melting.hide()
            self.mol_img.setPixmap(QPixmap(img))

    def smiles_action(self):
        if self.smiles_entry.isHidden():
            self.smiles_entry.show()
        else:
            self.smiles_entry.hide()


class EnterSmiles(QWidget):
    def __init__(self, mols, smiles, load):
        super().__init__()

        self.mols = mols
        self.smiles = smiles
        self.load = load

        self.resize(400, 100)

        # Set the overall layout
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.entry = QLineEdit()

        self.entry.setPlaceholderText('Enter SMILES string separated by commas')
        self.entry.returnPressed.connect(self.string_returned)

        self.layout.addWidget(self.entry)

    @pyqtSlot()
    def string_returned(self):
        self.mols.clear()
        self.smiles.clear()
        for mol in self.entry.text().replace(' ', '').split(','):
            try:
                self.mols.append(Chem.MolFromSmiles(mol))
                self.smiles.append(mol)
            except Exception as e:
                print(e)

        self.load()


app = QApplication(sys.argv)

window = MainWindow()
window.show()
sys.exit(app.exec_())
