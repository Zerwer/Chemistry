# Widget that allows manual entry of SMILES strings
from PyQt5.QtWidgets import *
from rdkit import Chem


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
        self.enter = QPushButton()

        self.entry.setPlaceholderText('Enter SMILES string separated by commas')
        self.entry.returnPressed.connect(self.string_returned)

        self.enter.setText('Enter')
        self.enter.clicked.connect(self.string_returned)

        self.layout.addWidget(self.entry)
        self.layout.addWidget(self.enter)

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
