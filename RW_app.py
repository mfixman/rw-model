#!/usr/bin/env python

from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QFormLayout, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
                             QPushButton, QTabWidget, QTableWidget, QVBoxLayout, QWidget, QStyleFactory)


class PavlovianApp(QDialog):
    def __init__(self, parent=None):
        super(PavlovianApp, self).__init__(parent)

        self.adaptive_types = ['linear', 'exponential', 'mack', 'hall', 'macknhall', 'dualV', 'lepelley', 'dualmack', 'hybrid']
        self.current_adaptive_type = None

        self.originalPalette = QApplication.palette()

        self.initUI()

        # Use a timer to ensure all widgets are updated after the window is shown
        QTimer.singleShot(100, self.updateWidgets)

    def initUI(self):
        styleComboBox = QComboBox()
        styleComboBox.addItems(QStyleFactory.keys())

        styleLabel = QLabel("&Style:")
        styleLabel.setBuddy(styleComboBox)

        self.useStylePaletteCheckBox = QCheckBox("&Use style's standard palette")
        self.useStylePaletteCheckBox.setChecked(True)

        disableWidgetsCheckBox = QCheckBox("&Disable widgets")

        self.createAdaptiveTypeGroupBox()
        self.createParametersGroupBox()
        self.createTableWidget()

        disableWidgetsCheckBox.toggled.connect(self.adaptiveTypeGroupBox.setDisabled)
        disableWidgetsCheckBox.toggled.connect(self.parametersGroupBox.setDisabled)
        disableWidgetsCheckBox.toggled.connect(self.tableTabWidget.setDisabled)

        mainLayout = QGridLayout()
        mainLayout.addWidget(self.tableTabWidget, 0, 0, 1, 2)  # Expanded to cover two columns
        mainLayout.addWidget(self.parametersGroupBox, 1, 0, 1, 1)  # Swapped position with adaptiveTypeGroupBox
        mainLayout.addWidget(self.adaptiveTypeGroupBox, 1, 1, 1, 1)  # Swapped position with parametersGroupBox
        self.setLayout(mainLayout)

        self.setWindowTitle("Pavlovian App")
        self.changeStyle('Windows')

    def changeStyle(self, styleName):
        QApplication.setStyle(QStyleFactory.create(styleName))
        self.changePalette()

    def changePalette(self):
        if self.useStylePaletteCheckBox.isChecked():
            QApplication.setPalette(QApplication.style().standardPalette())
        else:
            QApplication.setPalette(self.originalPalette)

    def createAdaptiveTypeGroupBox(self):
        self.adaptiveTypeGroupBox = QGroupBox("Adaptive Type")

        self.adaptivetypeComboBox = QComboBox(self)
        self.adaptivetypeComboBox.addItems(self.adaptive_types)
        self.adaptivetypeComboBox.activated.connect(self.changeAdaptiveType)

        layout = QVBoxLayout()
        layout.addWidget(self.adaptivetypeComboBox)
        layout.addStretch(1)
        self.adaptiveTypeGroupBox.setLayout(layout)

    def changeAdaptiveType(self):
        self.current_adaptive_type = self.adaptivetypeComboBox.currentText()
        print(self.current_adaptive_type)

        # Enable specific widgets based on the selected adaptive type
        widgets_to_enable = {
            'linear': ['alpha', 'lamda', 'beta'],
            'exponential': ['alpha', 'lamda', 'beta'],
            'mack': ['alpha', 'lamda', 'beta', 'thetaE', 'thetaI'],
            'hall': ['lamda', 'lamda', 'beta', 'gamma', 'thetaE', 'thetaI'],
            'macknhall': ['alpha', 'lamda', 'beta', 'gamma', 'window_size'],
            'dualV': ['alpha', 'lamda', 'beta', 'betan', 'gamma'],
            'lepelley': ['alpha', 'lamda', 'beta', 'betan', 'gamma', 'thetaE', 'thetaI'],
            'dualmack': ['alpha', 'lamda', 'beta', 'betan'],
            'hybrid': ['alpha', 'lamda', 'beta', 'betan', 'gamma', 'thetaE', 'thetaI'],
        }

        # Disable all widgets initially
        for key in ['alpha', 'lamda', 'beta', 'betan', 'gamma', 'thetaE', 'thetaI', 'window_size']:
            widget = getattr(self, f'{key}_box')
            widget.setDisabled(True)

        # Enable the widgets for the current adaptive type
        for key in widgets_to_enable[self.current_adaptive_type]:
            widget = getattr(self, f'{key}_box')
            widget.setDisabled(False)

    def createParametersGroupBox(self):
        self.parametersGroupBox = QGroupBox("Parameters")

        params = QFormLayout()
        self.alpha_Label = QLabel("Initial Alpha")
        self.alpha_box = QLineEdit()
        self.lamda_Label = QLabel("Lambda")
        self.lamda_box = QLineEdit()
        self.beta_Label = QLabel("Beta")
        self.beta_box = QLineEdit()
        self.betan_Label = QLabel("Inhibitory Beta")
        self.betan_box = QLineEdit()
        self.gamma_Label = QLabel("Gamma")
        self.gamma_box = QLineEdit()
        self.thetaE_Label = QLabel("Theta E")
        self.thetaE_box = QLineEdit()
        self.thetaI_Label = QLabel("Theta I")
        self.thetaI_box = QLineEdit()
        self.window_size_Label = QLabel("Window Size")
        self.window_size_box = QLineEdit()

        params.addRow(self.alpha_Label, self.alpha_box)
        params.addRow(self.lamda_Label, self.lamda_box)
        params.addRow(self.beta_Label, self.beta_box)
        params.addRow(self.betan_Label, self.betan_box)
        params.addRow(self.gamma_Label, self.gamma_box)
        params.addRow(self.thetaE_Label, self.thetaE_box)
        params.addRow(self.thetaI_Label, self.thetaI_box)
        params.addRow(self.window_size_Label, self.window_size_box)

        self.parametersGroupBox.setLayout(params)

    def createTableWidget(self):
        self.tableTabWidget = QTabWidget()

        tab1 = QWidget()
        self.tableWidget = QTableWidget(2, 1)  # 1 column, 2 rows

        addColumnButton = QPushButton("Add Column")
        addColumnButton.clicked.connect(self.addColumn)

        removeColumnButton = QPushButton("Remove Column")
        removeColumnButton.clicked.connect(self.removeColumn)

        addRowButton = QPushButton("Add Row")
        addRowButton.clicked.connect(self.addRow)

        removeRowButton = QPushButton("Remove Row")
        removeRowButton.clicked.connect(self.removeRow)

        buttonLayout = QVBoxLayout()
        buttonLayout.addWidget(addColumnButton)
        buttonLayout.addWidget(removeColumnButton)
        buttonLayout.addWidget(addRowButton)
        buttonLayout.addWidget(removeRowButton)
        buttonLayout.addStretch(1)

        tab1Layout = QHBoxLayout()
        tab1Layout.setContentsMargins(5, 5, 5, 5)
        tab1Layout.addWidget(self.tableWidget)
        tab1Layout.addLayout(buttonLayout)
        tab1.setLayout(tab1Layout)

        self.tableTabWidget.addTab(tab1, "&Table")

    def addColumn(self):
        currentColumnCount = self.tableWidget.columnCount()
        self.tableWidget.setColumnCount(currentColumnCount + 1)

    def removeColumn(self):
        currentColumnCount = self.tableWidget.columnCount()
        if currentColumnCount > 1:
            self.tableWidget.setColumnCount(currentColumnCount - 1)

    def addRow(self):
        currentRowCount = self.tableWidget.rowCount()
        self.tableWidget.setRowCount(currentRowCount + 1)

    def removeRow(self):
        currentRowCount = self.tableWidget.rowCount()
        if currentRowCount > 1:
            self.tableWidget.setRowCount(currentRowCount - 1)

    def updateWidgets(self):
        # Explicitly update and repaint widgets to ensure they are rendered correctly
        self.tableWidget.update()
        self.tableWidget.repaint()
        self.update()
        self.repaint()


if __name__ == '__main__':
    import sys

    app = QApplication(sys.argv)
    gallery = PavlovianApp()
    gallery.show()
    sys.exit(app.exec())
