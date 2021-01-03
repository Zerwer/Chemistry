from PyQt5.QtWidgets import *


class Graph:
    def __init__(self, window):
        self.window = window
        self.action = self.action()
        pass

    @property
    def action(self):
        raise NotImplementedError

    @action.setter
    def action(self, value):
        self._action = value

    def generate(self):
        pass


class Line(Graph):
    def action(self):
        act = QAction('Line', self.window)
        act.setStatusTip('Line Graph')
        return act

    def __init__(self, window):
        super().__init__(window)
        pass


class Scatter(Graph):
    def action(self):
        act = QAction('Scatter', self.window)
        act.setStatusTip('Scatter Plot')
        return act

    def __init__(self, window):
        super().__init__(window)
        pass


class AllGraphs:
    def __init__(self, window):
        self.line_graph = Line(window)
        self.scatter_graph = Scatter(window)

        self.actions = [self.line_graph.action,
                        self.scatter_graph.action]


class GraphOptions(QWidget):
    graph_types = ['line', 'scatter']
    descriptors = ['logP', 'solubility', 'melting', 'pKa']

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

        self.graph_name.setPlaceholderText('Enter name of graph')
        self.graph_type.addItems(self.graph_types)
        self.x_label.setPlaceholderText('Enter X-axis label')
        self.x_descriptor.addItems(self.descriptors)
        self.y_label.setPlaceholderText('Enter Y-axis label')
        self.y_descriptor.addItems(self.descriptors)
        self.enter.setText('Generate')
        self.enter.clicked.connect(self.generate_graph)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.layout.addWidget(self.graph_name)
        self.layout.addWidget(self.graph_type)
        self.layout.addWidget(self.x_label)
        self.layout.addWidget(self.x_descriptor)
        self.layout.addWidget(self.y_label)
        self.layout.addWidget(self.y_descriptor)
        self.layout.addWidget(self.enter)

    def generate_graph(self):
        pass
