# GUI app to interact with the models easily
import sys
from widgets.smiles import EnterSmiles
from widgets.properties import Properties
from widgets.graph import AllGraphs
from PyQt5.QtWidgets import *
from math import floor
from chemical_models import *

w, h = 800, 600


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.models = AllModels('logP', 'logS_logP', 'water_solubility',
                                'combined_solubility', 'melting_gse', 'pKa')

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
        self.properties = Properties(self.models, w, h)

        self.create_toolbar()
        self.create_search_bar()

        self.layout.addWidget(self.search_bar, 0, 0)
        self.layout.addWidget(self.table, 1, 0)
        self.layout.addWidget(self.properties, 1, 1)

        self.setCentralWidget(self.central_widget)

        self.properties.change_mol(None)
        self.mols_loaded()

    def create_search_bar(self):
        self.search_bar.setFixedWidth(floor(w / 2))

        self.search_bar.setPlaceholderText('Search...')
        self.search_bar.returnPressed.connect(self.string_searched)

    def create_toolbar(self):
        menu_bar = self.menuBar()
        menu_bar.setNativeMenuBar(True)

        graphs = AllGraphs(self)

        file_menu = menu_bar.addMenu('File')
        view_menu = menu_bar.addMenu('View')
        graph_menu = view_menu.addMenu('Graph')

        smiles_act = QAction('Enter SMILES', self)
        smiles_act.setShortcut('Ctrl+J')
        smiles_act.setStatusTip('Manually enter SMILES strings')

        smiles_act.triggered.connect(self.smiles_action)

        graph_menu.addActions(graphs.actions)

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
        self.properties.change_mol(self.displayed_mols[row])

    def smiles_action(self):
        if self.smiles_entry.isHidden():
            self.smiles_entry.show()
        else:
            self.smiles_entry.hide()


app = QApplication(sys.argv)

window = MainWindow()
window.show()
sys.exit(app.exec_())
