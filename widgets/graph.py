from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
import matplotlib.pyplot as plt
import io
import PIL
from PIL.ImageQt import ImageQt


# Data structure to combine all graphs
class AllGraphs:
    def __init__(self, window, models, mols):
        self.line_graph = Line(window, models, mols)
        self.scatter_graph = Scatter(window, models, mols)

        self.actions = [self.line_graph.action,
                        self.scatter_graph.action]


# All graphs inherent this general functionality
class Graph:
    def __init__(self, window, models, mols):
        self.window = window
        self.models = models
        self.mols = mols
        self.action = self.action()
        self.widget = self.Widget(self.generate_graph, self.models.descriptors.keys())

    # Widget displayed when selecting options for graph
    class Widget(QWidget):
        def __init__(self, generate_graph, descriptors):
            super().__init__()
            self.generate_graph = generate_graph
            self.descriptors = descriptors
            self.setWindowFlag(Qt.WindowStaysOnTopHint)

        # Allows graph specific arguments to be passed down
        def _generate_graph(self):
            pass

    @property
    def action(self):
        raise NotImplementedError

    @action.setter
    def action(self, value):
        self._action = value

    def generate_graph(self, **kwargs):
        pass

    def generate(self):
        if self.widget.isHidden():
            self.widget.show()
        else:
            self.widget.hide()


class DisplayGraph(QWidget):
    def __init__(self, title, img):
        super().__init__()
        self.setWindowTitle(title)

        self.resize(700, 500)

        self.graph = QLabel()
        self.save = QPushButton()

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.graph.setPixmap(img)
        self.save.setText('Save')

        layout.addWidget(self.graph)
        layout.addWidget(self.save)



class Line(Graph):
    def __init__(self, _win, _mod, _mol):
        super().__init__(_win, _mod, _mol)

    def action(self):
        act = QAction('Line', self.window)
        act.setStatusTip('Line Graph')
        act.triggered.connect(lambda: self.generate())
        return act

    def generate_graph(self):
        pass

    class Widget(Graph.Widget):
        def __init__(self, _, _d):
            super().__init__(_, _d)
            self.setWindowTitle('Select Options for Graph')

            self.resize(300, 200)

        def _generate_graph(self):
            pass


class Scatter(Graph):
    def __init__(self, _win, _mod, _mol):
        super().__init__(_win, _mod, _mol)

    def action(self):
        act = QAction('Scatter', self.window)
        act.setStatusTip('Scatter Plot')
        act.triggered.connect(lambda: self.generate())
        return act

    def generate_graph(self, graph_name, x_label,
                       x_descriptor, y_label, y_descriptor):
        x_abr_descriptor = self.models.descriptors[x_descriptor]
        y_abr_descriptor = self.models.descriptors[y_descriptor]
        x, y = [], []
        for mol in self.mols:
            properties = self.models.predict(mol, [x_abr_descriptor,
                                                   y_abr_descriptor])
            x.append(properties[x_abr_descriptor])
            y.append(properties[y_abr_descriptor])

        plt.scatter(x, y, s=1)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(graph_name)

        buf = io.BytesIO()
        plt.savefig(buf)
        plt.close()
        buf.seek(0)
        img = ImageQt(PIL.Image.open(buf))
        self.widget.close()
        self.window.display_graph(graph_name, QPixmap.fromImage(img))

    class Widget(Graph.Widget):
        def __init__(self, _, _d):
            super().__init__(_, _d)
            self.setWindowTitle('Select Options for Graph')

            self.resize(300, 200)

            self.graph_name = QLineEdit()
            self.x_label = QLineEdit()
            self.x_descriptor = QComboBox()
            self.y_label = QLineEdit()
            self.y_descriptor = QComboBox()
            self.enter = QPushButton()

            self.graph_name.setPlaceholderText('Enter name of graph')
            self.x_label.setPlaceholderText('Enter X-axis label')
            self.x_descriptor.addItems(self.descriptors)
            self.y_label.setPlaceholderText('Enter Y-axis label')
            self.y_descriptor.addItems(self.descriptors)
            self.enter.setText('Generate')
            self.enter.clicked.connect(self._generate_graph)

            layout = QVBoxLayout()
            self.setLayout(layout)

            layout.addWidget(self.graph_name)
            layout.addWidget(self.x_label)
            layout.addWidget(self.x_descriptor)
            layout.addWidget(self.y_label)
            layout.addWidget(self.y_descriptor)
            layout.addWidget(self.enter)

        def _generate_graph(self):
            kwargs = {'graph_name': self.graph_name.text(),
                      'x_label': self.x_label.text(),
                      'x_descriptor': str(self.x_descriptor.currentText()),
                      'y_label': self.y_label.text(),
                      'y_descriptor': str(self.y_descriptor.currentText())}
            self.generate_graph(**kwargs)
