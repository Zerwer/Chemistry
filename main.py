"""
GUI app to interact with the models easily
"""
import sys
from PIL.ImageQt import ImageQt
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QLabel, QWidget, QVBoxLayout, QLineEdit
from PyQt5.QtGui import QImage, QPixmap, QValidator, QDoubleValidator
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors
from chemical_models import AtomPairSolubility, LogP, LogPSolubility, CombinedSolubility

# Load models
logP_model = LogP('logP')


class MainWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.mol_img = QLabel()
        self.logP = QLabel()

        self.smiles = QLineEdit()
        self.smiles.returnPressed.connect(self.string_returned)

        self.layout = QVBoxLayout()
        self.layout.addWidget(self.mol_img)
        self.layout.addWidget(self.logP)
        self.layout.addWidget(self.smiles)
        self.setLayout(self.layout)

    @pyqtSlot()
    def string_returned(self):
        try:
            mol = Chem.MolFromSmiles(self.smiles.text())
            img = QImage(ImageQt(Draw.MolToImage(mol)))
            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
            logP = logP_model.run(fp)

            self.logP.setText('LogP: ' + str(logP))
            self.mol_img.setPixmap(QPixmap(img))
        except Exception as e:
            print(e)


app = QApplication(sys.argv)

widget = MainWidget()
widget.resize(800, 600)
widget.show()
sys.exit(app.exec_())
