# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow.ui'
#
# Created: Mon Sep  8 10:58:47 2014
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.pyplot import gcf, setp

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(958, 1255)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(20, 10, 101, 41))
        self.label.setObjectName(_fromUtf8("label"))
        self.pushButton = QtGui.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(560, 20, 98, 27))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.graphicsView = QtGui.QGraphicsView(self.centralwidget)
        self.graphicsView.setGeometry(QtCore.QRect(40, 230, 721, 411))
        self.graphicsView.setObjectName(_fromUtf8("graphicsView"))
        self.pushButton_2 = QtGui.QPushButton(self.centralwidget)
        self.pushButton_2.setGeometry(QtCore.QRect(560, 70, 98, 27))
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(20, 70, 66, 17))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.lineEdit = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit.setGeometry(QtCore.QRect(130, 70, 381, 27))
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.lineEdit_2 = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit_2.setGeometry(QtCore.QRect(130, 20, 381, 27))
        self.lineEdit_2.setObjectName(_fromUtf8("lineEdit_2"))
        self.horizontalSlider = QtGui.QSlider(self.centralwidget)
        self.horizontalSlider.setGeometry(QtCore.QRect(40, 650, 251, 29))
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setTickPosition(QtGui.QSlider.TicksAbove)
        self.horizontalSlider.setObjectName(_fromUtf8("horizontalSlider"))
        self.verticalSlider = QtGui.QSlider(self.centralwidget)
        self.verticalSlider.setGeometry(QtCore.QRect(780, 230, 29, 121))
        self.verticalSlider.setOrientation(QtCore.Qt.Vertical)
        self.verticalSlider.setTickPosition(QtGui.QSlider.TicksAbove)
        self.verticalSlider.setObjectName(_fromUtf8("verticalSlider"))
        self.comboBox = QtGui.QComboBox(self.centralwidget)
        self.comboBox.setGeometry(QtCore.QRect(20, 130, 491, 27))
        self.comboBox.setObjectName(_fromUtf8("comboBox"))
        self.pushButton_3 = QtGui.QPushButton(self.centralwidget)
        self.pushButton_3.setGeometry(QtCore.QRect(560, 130, 98, 27))
        self.pushButton_3.setObjectName(_fromUtf8("pushButton_3"))
        self.lineEdit_3 = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit_3.setGeometry(QtCore.QRect(140, 680, 151, 27))
        self.lineEdit_3.setObjectName(_fromUtf8("lineEdit_3"))
        self.lineEdit_4 = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit_4.setGeometry(QtCore.QRect(602, 680, 161, 27))
        self.lineEdit_4.setObjectName(_fromUtf8("lineEdit_4"))
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(40, 680, 101, 17))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(490, 680, 101, 17))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.filterTab = QtGui.QTabWidget(self.centralwidget)
        self.filterTab.setGeometry(QtCore.QRect(40, 800, 731, 351))
        self.filterTab.setObjectName(_fromUtf8("filterTab"))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.comboBox_2 = QtGui.QComboBox(self.tab)
        self.comboBox_2.setGeometry(QtCore.QRect(150, 40, 421, 27))
        self.comboBox_2.setObjectName(_fromUtf8("comboBox_2"))
        self.label_7 = QtGui.QLabel(self.tab)
        self.label_7.setGeometry(QtCore.QRect(20, 40, 81, 17))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.label_8 = QtGui.QLabel(self.tab)
        self.label_8.setGeometry(QtCore.QRect(40, 90, 111, 16))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.label_9 = QtGui.QLabel(self.tab)
        self.label_9.setGeometry(QtCore.QRect(390, 90, 111, 17))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.label_10 = QtGui.QLabel(self.tab)
        self.label_10.setGeometry(QtCore.QRect(50, 140, 101, 16))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.label_11 = QtGui.QLabel(self.tab)
        self.label_11.setGeometry(QtCore.QRect(350, 140, 161, 17))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.pushButton_4 = QtGui.QPushButton(self.tab)
        self.pushButton_4.setGeometry(QtCore.QRect(600, 40, 98, 27))
        self.pushButton_4.setObjectName(_fromUtf8("pushButton_4"))
        self.label_19 = QtGui.QLabel(self.tab)
        self.label_19.setGeometry(QtCore.QRect(10, 190, 141, 16))
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.label_20 = QtGui.QLabel(self.tab)
        self.label_20.setGeometry(QtCore.QRect(40, 240, 111, 17))
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.label_21 = QtGui.QLabel(self.tab)
        self.label_21.setGeometry(QtCore.QRect(410, 190, 101, 17))
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.lineEdit_5 = QtGui.QLineEdit(self.tab)
        self.lineEdit_5.setGeometry(QtCore.QRect(160, 90, 181, 27))
        self.lineEdit_5.setObjectName(_fromUtf8("lineEdit_5"))
        self.lineEdit_6 = QtGui.QLineEdit(self.tab)
        self.lineEdit_6.setGeometry(QtCore.QRect(530, 90, 181, 27))
        self.lineEdit_6.setObjectName(_fromUtf8("lineEdit_6"))
        self.lineEdit_7 = QtGui.QLineEdit(self.tab)
        self.lineEdit_7.setGeometry(QtCore.QRect(160, 140, 181, 27))
        self.lineEdit_7.setObjectName(_fromUtf8("lineEdit_7"))
        self.lineEdit_8 = QtGui.QLineEdit(self.tab)
        self.lineEdit_8.setGeometry(QtCore.QRect(160, 190, 181, 27))
        self.lineEdit_8.setObjectName(_fromUtf8("lineEdit_8"))
        self.lineEdit_9 = QtGui.QLineEdit(self.tab)
        self.lineEdit_9.setGeometry(QtCore.QRect(160, 240, 181, 27))
        self.lineEdit_9.setObjectName(_fromUtf8("lineEdit_9"))
        self.comboBox_4 = QtGui.QComboBox(self.tab)
        self.comboBox_4.setGeometry(QtCore.QRect(530, 140, 181, 27))
        self.comboBox_4.setObjectName(_fromUtf8("comboBox_4"))
        self.comboBox_5 = QtGui.QComboBox(self.tab)
        self.comboBox_5.setGeometry(QtCore.QRect(530, 190, 181, 27))
        self.comboBox_5.setObjectName(_fromUtf8("comboBox_5"))
        self.checkBox = QtGui.QCheckBox(self.tab)
        self.checkBox.setGeometry(QtCore.QRect(420, 240, 97, 22))
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.pushButton_6 = QtGui.QPushButton(self.tab)
        self.pushButton_6.setGeometry(QtCore.QRect(590, 250, 121, 51))
        self.pushButton_6.setObjectName(_fromUtf8("pushButton_6"))
        self.filterTab.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.lineEdit_timeInterval = QtGui.QLineEdit(self.tab_2)
        self.lineEdit_timeInterval.setGeometry(QtCore.QRect(180, 30, 191, 27))
        self.lineEdit_timeInterval.setObjectName(_fromUtf8("lineEdit_timeInterval"))
        self.label_6 = QtGui.QLabel(self.tab_2)
        self.label_6.setGeometry(QtCore.QRect(20, 30, 101, 17))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.pushButton_5 = QtGui.QPushButton(self.tab_2)
        self.pushButton_5.setGeometry(QtCore.QRect(500, 30, 121, 51))
        self.pushButton_5.setObjectName(_fromUtf8("pushButton_5"))
        self.filterTab.addTab(self.tab_2, _fromUtf8(""))
        self.verticalSlider_2 = QtGui.QSlider(self.centralwidget)
        self.verticalSlider_2.setGeometry(QtCore.QRect(780, 520, 29, 121))
        self.verticalSlider_2.setOrientation(QtCore.Qt.Vertical)
        self.verticalSlider_2.setTickPosition(QtGui.QSlider.TicksBelow)
        self.verticalSlider_2.setObjectName(_fromUtf8("verticalSlider_2"))
        self.horizontalSlider_2 = QtGui.QSlider(self.centralwidget)
        self.horizontalSlider_2.setGeometry(QtCore.QRect(510, 650, 251, 29))
        self.horizontalSlider_2.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_2.setTickPosition(QtGui.QSlider.TicksAbove)
        self.horizontalSlider_2.setObjectName(_fromUtf8("horizontalSlider_2"))
        self.lineEdit_timeInterval_2 = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit_timeInterval_2.setGeometry(QtCore.QRect(110, 750, 261, 27))
        self.lineEdit_timeInterval_2.setObjectName(_fromUtf8("lineEdit_timeInterval_2"))
        self.label_12 = QtGui.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(50, 750, 141, 17))
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.lineEdit_timeInterval_3 = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit_timeInterval_3.setGeometry(QtCore.QRect(560, 750, 211, 27))
        self.lineEdit_timeInterval_3.setObjectName(_fromUtf8("lineEdit_timeInterval_3"))
        self.label_13 = QtGui.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(420, 750, 141, 17))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 958, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.figure = Figure((4.0, 3.0), dpi=100)
        self.theplot = self.figure.add_subplot(111)
        self.graphicsView = FigureCanvas(self.figure)
        self.graphicsView.setParent(self.centralwidget)
        self.graphicsView.setGeometry(QtCore.QRect(40, 230, 721, 411))
        self.graphicsView.setObjectName(_fromUtf8("graphicsView"))

        self.retranslateUi(MainWindow)
        self.filterTab.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Run Number", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("MainWindow", "Load", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_2.setText(QtGui.QApplication.translate("MainWindow", "Load", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "File Name", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_3.setText(QtGui.QApplication.translate("MainWindow", "Use", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Starting Time", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "Stopping Time", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("MainWindow", "Sample Log", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("MainWindow", "Minimum Value", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("MainWindow", "Maximum Value", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("MainWindow", "Log Step Size", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("MainWindow", "Value Change Direction", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_4.setText(QtGui.QApplication.translate("MainWindow", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.label_19.setText(QtGui.QApplication.translate("MainWindow", "Log Value Tolerance", None, QtGui.QApplication.UnicodeUTF8))
        self.label_20.setText(QtGui.QApplication.translate("MainWindow", "Time Tolerance", None, QtGui.QApplication.UnicodeUTF8))
        self.label_21.setText(QtGui.QApplication.translate("MainWindow", "Log Boundary", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBox.setText(QtGui.QApplication.translate("MainWindow", "Fast Log", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_6.setText(QtGui.QApplication.translate("MainWindow", "Filter", None, QtGui.QApplication.UnicodeUTF8))
        self.filterTab.setTabText(self.filterTab.indexOf(self.tab), QtGui.QApplication.translate("MainWindow", "Filter By Log", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("MainWindow", "Time Interval", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_5.setText(QtGui.QApplication.translate("MainWindow", "Filter", None, QtGui.QApplication.UnicodeUTF8))
        self.filterTab.setTabText(self.filterTab.indexOf(self.tab_2), QtGui.QApplication.translate("MainWindow", "Filter By Time", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("MainWindow", "Title", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("MainWindow", "Output Workspace", None, QtGui.QApplication.UnicodeUTF8))

