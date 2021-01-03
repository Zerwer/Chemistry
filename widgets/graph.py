from PyQt5.QtWidgets import *


class GraphOptions(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Select Options for Graph')

        self.resize(300, 200)

        self.graph_name = QLineEdit()
        self.graph_type = QComboBox()
        self.x_label = QLineEdit()
        self.x_descriptor = QComboBox()
        self.y_label = QLineEdit()
        self.y_descriptor = QComboBox()
        self.enter = QPushButton()

        self.enter.setText('Enter')
        # self.enter.clicked.connect()

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.layout.addWidget(self.graph_name)
        self.layout.addWidget(self.graph_type)
        self.layout.addWidget(self.x_label)
        self.layout.addWidget(self.x_descriptor)
        self.layout.addWidget(self.y_label)
        self.layout.addWidget(self.y_descriptor)
        self.layout.addWidget(self.enter)
