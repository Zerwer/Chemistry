from PyQt5.QtWidgets import *


class GraphOptions(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Select Options for Graph')

        self.resize(200, 200)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)