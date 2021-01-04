from PIL.ImageQt import ImageQt
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QImage, QPixmap
from math import floor
from rdkit.Chem import Draw


class Properties(QWidget):
    def __init__(self, models, w, h):
        super().__init__()

        self.w, self.h = w, h

        self.models = models

        self.mol = None

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.mol_img = QLabel()
        self.logP = QLabel()
        self.solubility = QLabel()
        self.melting = QLabel()
        self.pKa = QLabel()

        self.layout.addWidget(self.mol_img)
        self.layout.addWidget(self.logP)
        self.layout.addWidget(self.solubility)
        self.layout.addWidget(self.melting)
        self.layout.addWidget(self.pKa)

    def change_mol(self, mol):
        self.mol = mol

        if self.mol is None:
            self.logP.hide()
            self.solubility.hide()
            self.melting.hide()
            self.pKa.hide()

            self.mol_img.setText('No Molecule Selected!')
        else:
            img = QImage(ImageQt(Draw.MolToImage(self.mol, size=(floor(self.w/2),
                                                                 floor(self.h/2)))))

            predictions = self.models.predict(self.mol, ['logP', 'sol',
                                                         'mp', 'pka'])

            self.logP.show()
            self.solubility.show()
            self.melting.show()
            self.pKa.show()

            self.logP.setText('LogP: ' + str(round(predictions['logP'], 2)))
            self.solubility.setText('Water Solubility(mg/mL): ' +
                                    str(round(predictions['sol'], 2)))
            self.melting.setText('Melting Point(C): ' + str(round(predictions['mp'], 2)))
            self.pKa.setText('pKa: ' + str(round(predictions['pka'], 2)))

            self.mol_img.setPixmap(QPixmap(img))
