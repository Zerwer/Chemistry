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


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Stores both the RDKit Mol and SMILES
        self.mols = []
        self.smiles = []

        self.toolbar = QToolBar()
        self.widget = MainWidget()

        self.smiles = EnterSmiles(self.mols, self.smiles)

        self.setCentralWidget(self.widget)
        self.create_toolbar()

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
        if len(self.mols) == 1:
            self.calculate_properties(self.mols[0], self.smiles[0])
        elif len(self.mols) <= 35:
            fps = []
            for mol in self.mols:
                fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
                fps.append((fingerprint, mol))
            dim = int(ceil(sqrt(len(self.mols))))

            img = QImage(ImageQt(mol_similarity_grid(fps, (int(700/dim),
                                                     int(700/dim)), dim)))
            self.mol_img.setPixmap(QPixmap(img))

    def smiles_action(self):
        if self.smiles.isHidden():
            self.smiles.show()
        else:
            self.smiles.hide()


class EnterSmiles(QWidget):
    def __init__(self, mols, smiles):
        super().__init__()

        self.mols = mols
        self.smiles = smiles

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
        self.mols = []
        for mol in self.entry.text().replace(' ', '').split(','):
            try:
                self.mols.append(Chem.MolFromSmiles(mol))
                self.smiles.append(mol)
            except Exception as e:
                print(e)

        self.mols_loaded()


# Main widget of the application
class MainWidget(QTabWidget):
    def __init__(self):
        QTabWidget.__init__(self)

        # Stores both the RDKit Mol and SMILES
        self.mols = []
        self.smiles = []

        # Set the overall layout
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        # Initiate all elements here
        self.selection_instruction = QLabel()
        self.select_smiles_file = QPushButton('Select File')
        self.mol_img = QLabel()
        self.mol_name = QLabel()
        self.logP = QLabel()
        self.solubility = QLabel()
        self.melting = QLabel()
        self.entry = QLineEdit()

        # Configure widgets
        self.entry.setPlaceholderText('Enter SMILES string separated by commas')
        self.setWindowTitle('Chemistry')

        # Connect widgets to slots
        self.select_smiles_file.pressed.connect(self.select_file)
        self.entry.returnPressed.connect(self.string_returned)

        # Create and add tabs here
        self.mol_selection = QWidget()
        self.mol_display = QWidget()
        self.mol_graph = QWidget()
        self.mol_properties = QWidget()

        self.addTab(self.mol_selection, 'Select')
        self.addTab(self.mol_display, 'View')
        self.addTab(self.mol_graph, 'Graph')
        self.addTab(self.mol_properties, 'Properties')

        # Initiate tabs
        self.mol_selection_ui()
        self.mol_display_ui()
        self.mol_properties_ui()

    # UI layout functions got here:
    def mol_selection_ui(self):
        layout = QVBoxLayout()
        layout.addWidget(self.entry)
        layout.addWidget(self.selection_instruction)
        layout.addWidget(self.select_smiles_file)
        self.selection_instruction.setText('Or select a file containing a list '
                                           'of SMILES strings')
        self.mol_selection.setLayout(layout)

    def mol_display_ui(self):
        layout = QVBoxLayout()
        layout.addWidget(self.mol_img)
        self.mol_display.setLayout(layout)

    def mol_properties_ui(self):
        layout = QVBoxLayout()
        layout.addWidget(self.mol_name)
        layout.addWidget(self.logP)
        layout.addWidget(self.solubility)
        layout.addWidget(self.melting)
        self.mol_properties.setLayout(layout)

    # Other functions:

    def calculate_properties(self, mol, smile):
        img = QImage(ImageQt(Draw.MolToImage(mol, size=(700, 700))))
        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
        logP = logP_model.run(fp)
        logP_sol = logP_solubility_model.run(logP)
        atom_pair_sol = atom_pair_sol_model.run(fp)
        combined_sol = combined_model.run(mol, logP,
                                          logP_sol, atom_pair_sol)

        mg_ml_sol = logs_to_mg_ml(combined_sol, mol)
        mp = melting_point_model.run(combined_sol, logP)

        self.mol_name.setText(smile)
        self.logP.setText('LogP: ' + str(round(logP, 2)))
        self.solubility.setText('Water Solubility(mg/mL): ' +
                                str(round(mg_ml_sol, 2)))
        self.melting.setText('Melting Point(C): ' + str(round(mp, 2)))
        self.mol_img.setPixmap(QPixmap(img))

    def mols_loaded(self):
        if len(self.mols) == 1:
            self.calculate_properties(self.mols[0], self.smiles[0])
        elif len(self.mols) <= 35:
            fps = []
            for mol in self.mols:
                fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
                fps.append((fingerprint, mol))
            dim = int(ceil(sqrt(len(self.mols))))

            img = QImage(ImageQt(mol_similarity_grid(fps, (int(700/dim),
                                                     int(700/dim)), dim)))
            self.mol_img.setPixmap(QPixmap(img))

    @pyqtSlot()
    def string_returned(self):
        self.mols = []
        for mol in self.entry.text().replace(' ', '').split(','):
            try:
                self.mols.append(Chem.MolFromSmiles(mol))
                self.smiles.append(mol)
            except Exception as e:
                print(e)

        self.mols_loaded()

    def select_file(self):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)

        if dlg.exec_():
            filenames = dlg.selectedFiles()
            file = open(filenames[0], 'r')

            # for line in file.readlines()
            # with f:
            #     data = f.read()
            #     self.contents.setText(data)


app = QApplication(sys.argv)

window = MainWindow()
window.resize(800, 600)
window.show()
sys.exit(app.exec_())
