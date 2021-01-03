from PIL.ImageQt import ImageQt
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QImage, QPixmap
from math import floor
from rdkit.Chem import Draw, MACCSkeys
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from chemical_models import *
from functions import logs_to_mg_ml


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

            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(self.mol)
            avalon = GetAvalonFP(mol)
            maacs = MACCSkeys.GenMACCSKeys(mol)

            logP = self.models.logP_model.run(fp)
            logP_sol = self.models.logP_solubility_model.run(logP)
            atom_pair_sol = self.models.atom_pair_sol_model.run(fp)
            combined_sol = self.models.combined_model.run(self.mol, logP,
                                              logP_sol, atom_pair_sol)
            mg_ml_sol = logs_to_mg_ml(combined_sol, self.mol)
            mp = self.models.melting_point_model.run(combined_sol, logP)
            pka = self.models.pKa_model.run(avalon + maacs + fp)

            self.logP.show()
            self.solubility.show()
            self.melting.show()
            self.pKa.show()

            self.logP.setText('LogP: ' + str(round(logP, 2)))
            self.solubility.setText('Water Solubility(mg/mL): ' +
                                    str(round(mg_ml_sol, 2)))
            self.melting.setText('Melting Point(C): ' + str(round(mp, 2)))
            self.pKa.setText('pKa: ' + str(round(pka, 2)))

            self.mol_img.setPixmap(QPixmap(img))
