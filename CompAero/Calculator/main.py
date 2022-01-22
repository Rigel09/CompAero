from CompAero.Calculator.CalculatorUI import Ui_MainWindow
from PyQt5.QtWidgets import QMainWindow, QWidget
from PyQt5 import QtCore, QtWidgets
from typing import Optional, Union
import sys


class UI(QMainWindow, Ui_MainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setupUi(self)
        self.show()


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    _ = UI()
    app.exec_()
