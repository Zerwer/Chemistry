from PyQt5.QtWidgets import *


# Data structure to combine all graphs
class AllGraphs:
    def __init__(self, window):
        self.line_graph = Line(window)
        self.scatter_graph = Scatter(window)

        self.actions = [self.line_graph.action,
                        self.scatter_graph.action]


# All graphs inherent this general functionality
class Graph:
    def __init__(self, window):
        self.window = window
        self.action = self.action()
        self.widget = self.Widget(self.generate_graph)

    class Widget(QWidget):
        descriptors = ['logP', 'solubility', 'melting', 'pKa']

        def __init__(self, generate_graph):
            super().__init__()
            self.generate_graph = generate_graph

    @property
    def action(self):
        raise NotImplementedError

    @action.setter
    def action(self, value):
        self._action = value

    def generate_graph(self):
        pass

    def generate(self):
        if self.widget.isHidden():
            self.widget.show()
        else:
            self.widget.hide()


class Line(Graph):
    def __init__(self, window):
        super().__init__(window)

    def action(self):
        act = QAction('Line', self.window)
        act.setStatusTip('Line Graph')
        act.triggered.connect(lambda: self.generate())
        return act

    class Widget(Graph.Widget):
        def __init__(self, _):
            super().__init__(_)
            self.setWindowTitle('Select Options for Graph')

            self.resize(300, 200)


class Scatter(Graph):
    def __init__(self, window):
        super().__init__(window)

    def action(self):
        act = QAction('Scatter', self.window)
        act.setStatusTip('Scatter Plot')
        act.triggered.connect(lambda: self.generate())
        return act

    class Widget(Graph.Widget):
        def __init__(self, _):
            super().__init__(_)
            self.setWindowTitle('Select Options for Graph')

            self.resize(300, 200)

            graph_name = QLineEdit()
            x_label = QLineEdit()
            x_descriptor = QComboBox()
            y_label = QLineEdit()
            y_descriptor = QComboBox()
            enter = QPushButton()

            graph_name.setPlaceholderText('Enter name of graph')
            x_label.setPlaceholderText('Enter X-axis label')
            x_descriptor.addItems(self.descriptors)
            y_label.setPlaceholderText('Enter Y-axis label')
            y_descriptor.addItems(self.descriptors)
            enter.setText('Generate')
            enter.clicked.connect(self.generate_graph)

            layout = QVBoxLayout()
            self.setLayout(layout)

            layout.addWidget(graph_name)
            layout.addWidget(x_label)
            layout.addWidget(x_descriptor)
            layout.addWidget(y_label)
            layout.addWidget(y_descriptor)
            layout.addWidget(enter)
