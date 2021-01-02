# GUI app to interact with the models easily
import sys
from PIL.ImageQt import ImageQt
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QImage, QPixmap
from math import floor
from rdkit.Chem import Draw
from chemical_models import *
from functions import logs_to_mg_ml

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
        self.layout = QGridLayout()
        self.central_widget.setLayout(self.layout)

        # Stores both the RDKit Mol and SMILES
        self.mols = []
        self.smiles = []
        # Mols that should be displayed
        self.displayed_mols = []
        self.displayed_smiles = []

        self.toolbar = QToolBar()
        self.search_bar = QLineEdit()
        self.table = QTableWidget()

        self.table.setFixedWidth(floor(w/2))

        self.smiles_entry = EnterSmiles(self.mols, self.smiles,
                                        self.displayed_mols,
                                        self.displayed_smiles, self.mols_loaded)
        self.properties = Properties()

        self.create_toolbar()
        self.create_search_bar()

        self.layout.addWidget(self.search_bar, 0, 0)
        self.layout.addWidget(self.table, 1, 0)
        self.layout.addWidget(self.properties, 1, 1)

        self.setCentralWidget(self.central_widget)

        self.properties.change_mol()
        self.mols_loaded()

    def create_search_bar(self):
        self.search_bar.setFixedWidth(floor(w / 2))

        self.search_bar.setPlaceholderText('Search...')
        self.search_bar.returnPressed.connect(self.string_searched)

    def create_toolbar(self):
        menu_bar = self.menuBar()
        menu_bar.setNativeMenuBar(True)

        file_menu = menu_bar.addMenu(' &File')

        smiles_act = QAction('Enter SMILES', self)

        smiles_act.setShortcut('Ctrl+J')
        smiles_act.setStatusTip('Manually enter SMILES strings')
        smiles_act.triggered.connect(self.smiles_action)

        file_menu.addAction(smiles_act)
        file_menu.addAction('Open...')

        self.setMenuBar(menu_bar)

    def mols_loaded(self):
        self.table.clear()

        if len(self.displayed_mols) == 0:
            return

        self.table.setRowCount(len(self.displayed_mols))
        self.table.setColumnCount(1)

        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setVisible(False)

        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)

        self.table.setColumnWidth(0, floor(w/2))

        self.table.cellClicked.connect(self.item_selected)
        self.table.cellActivated.connect(self.item_selected)
        self.table.cellEntered.connect(self.item_selected)

        for i, smile in enumerate(self.displayed_smiles):
            self.table.setRowHeight(i, 50)
            item = QTableWidgetItem(smile)
            self.table.setItem(i, 0, item)

    # Search checks for SMILES similarities
    def string_searched(self):
        self.displayed_mols.clear()
        self.displayed_smiles.clear()

        for mol, smile in zip(self.mols, self.smiles):
            if self.search_bar.text() in smile:
                self.displayed_mols.append(mol)
                self.displayed_smiles.append(smile)

        self.mols_loaded()

    def item_selected(self, row, _):
        self.properties.mol = self.displayed_mols[row]
        self.properties.change_mol()

    def smiles_action(self):
        if self.smiles_entry.isHidden():
            self.smiles_entry.show()
        else:
            self.smiles_entry.hide()


class Properties(QWidget):
    def __init__(self):
        super().__init__()

        self.mol = None

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.mol_img = QLabel()
        self.logP = QLabel()
        self.solubility = QLabel()
        self.melting = QLabel()

        self.layout.addWidget(self.mol_img)
        self.layout.addWidget(self.logP)
        self.layout.addWidget(self.solubility)
        self.layout.addWidget(self.melting)

    def change_mol(self):
        if self.mol is None:
            self.logP.hide()
            self.solubility.hide()
            self.melting.hide()

            self.mol_img.setText('No Molecule Selected!')
        else:
            img = QImage(ImageQt(Draw.MolToImage(self.mol, size=(floor(w/2), floor(h/2)))))

            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(self.mol)
            logP = logP_model.run(fp)
            logP_sol = logP_solubility_model.run(logP)
            atom_pair_sol = atom_pair_sol_model.run(fp)
            combined_sol = combined_model.run(self.mol, logP,
                                              logP_sol, atom_pair_sol)
            mg_ml_sol = logs_to_mg_ml(combined_sol, self.mol)
            mp = melting_point_model.run(combined_sol, logP)

            self.logP.show()
            self.solubility.show()
            self.melting.show()

            self.logP.setText('LogP: ' + str(round(logP, 2)))
            self.solubility.setText('Water Solubility(mg/mL): ' +
                                    str(round(mg_ml_sol, 2)))
            self.melting.setText('Melting Point(C): ' + str(round(mp, 2)))

            self.mol_img.setPixmap(QPixmap(img))


class EnterSmiles(QWidget):
    def __init__(self, mols, smiles, displayed_mols, displayed_smiles, load):
        super().__init__()

        self.setWindowTitle('Enter SMILES String')

        self.mols = mols
        self.smiles = smiles
        self.displayed_mols = displayed_mols
        self.displayed_smiles = displayed_smiles
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
        self.displayed_mols.clear()
        self.displayed_smiles.clear()

        for mol in self.entry.text().replace(' ', '').split(','):
            try:
                rd_mol = Chem.MolFromSmiles(mol)
                self.mols.append(rd_mol)
                self.displayed_mols.append(rd_mol)
                self.smiles.append(mol)
                self.displayed_smiles.append(mol)
            except Exception as e:
                print(e)

        self.load()
        self.close()


app = QApplication(sys.argv)

window = MainWindow()
window.show()
sys.exit(app.exec_())
