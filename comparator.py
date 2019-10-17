#!/usr/bin/env python3

import os
import sys
import random
import argparse
import numpy as np
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QLabel, QVBoxLayout, QHBoxLayout, QSizePolicy
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class App(QWidget):
    def __init__(self, files=None, edge_snip=False, baseline_adjust=False):
        super().__init__()

        self.edge_snip = edge_snip
        self.baseline_adjust = baseline_adjust

        if len(files) >= 1:
            self.file1 = files[0]
        else:
            self.file1 = None

        if len(files) >= 2:
            self.file2 = files[1]
        else:
            self.file2 = None

        self.title = 'Comparador'
        self.initUI()

        return

    def initUI(self):
        self.setWindowTitle(self.title)

        self.resize(600, 400)

        label1 = QLabel()
        openButton1 = QPushButton('File 1...')
        openButton1.clicked.connect(self.openFile1)
        self.label1 = label1

        label2 = QLabel()
        openButton2 = QPushButton('File 2...')
        openButton2.clicked.connect(self.openFile2)
        self.label2 = label2

        calcButton = QPushButton('Calculate similarity')
        calcButton.clicked.connect(self.calcSimilarity)
        resultLabel = QLabel('Select files to be compared')
        self.resultLabel = resultLabel

        self.updateFile1()
        self.updateFile2()


        sidePanel = QVBoxLayout()

        sidePanel.addWidget(label1)
        sidePanel.addWidget(openButton1)
        sidePanel.addWidget(label2)
        sidePanel.addWidget(openButton2)

        sidePanel.addStretch()

        sidePanel.addWidget(calcButton)
        sidePanel.addWidget(resultLabel)

        rootLayout = QHBoxLayout()

        rootLayout.addLayout(sidePanel)


        fig = plt.figure()
        self.fig = fig

        axes = fig.add_subplot(111)
        self.axes = axes

        axes.invert_xaxis()

        canvas = FigureCanvasToolbar(fig)
        self.canvas = canvas

        rootLayout.addWidget(canvas)

        self.setLayout(rootLayout)

        self.show()
        return

    @pyqtSlot()
    def openFile1(self):
        self.file1 = self.openFileDialog()
        self.updateFile1()
        return

    @pyqtSlot()
    def openFile2(self):
        self.file2 = self.openFileDialog()
        self.updateFile2()
        return

    @pyqtSlot()
    def calcSimilarity(self):
        sp1 = Spectrum.open(self.file1)
        sp2 = Spectrum.open(self.file2)

        #if False:
        if self.baseline_adjust:
            sp1 = sp1.baseline()
            sp2 = sp2.baseline()

        #if False:
        if self.edge_snip:
            sp1 = sp1.snip_edges()
            sp2 = sp2.snip_edges()

        sp1 = sp1.scaled2()
        sp2 = sp2.scaled2()

        result = Spectrum.similarity(sp1, sp2)

        print("aaa", result)

        self.resultLabel.setText('{:.5%}'.format(result))

        self.axes.clear()
        self.axes.invert_xaxis()
        self.axes.plot(sp1.x, sp1.real)
        self.axes.plot(sp2.x, sp2.real)
        self.canvas.draw()

        return

    def updateFile1(self):
        if self.file1 == None:
            self.label1.setText('No file selected')
        else:
            self.label1.setText(self.file1)
        return

    def updateFile2(self):
        if self.file2 == None:
            self.label2.setText('No file selected')
        else:
            self.label2.setText(self.file2)
        return

    def openFileDialog(self):
        file = QFileDialog.getOpenFileName(self, "Open Image", "", "CSV Files (*.asc)")

        if file[0] != '':
            return file[0]
        else:
            return None

class FigureCanvasToolbar(QWidget):
    def __init__(self, figure):
        super().__init__()

        layout = QVBoxLayout()
        self.setLayout(layout)

        canvas = FigureCanvas(figure)
        self.canvas = canvas

        toolbar = NavigationToolbar(canvas, QWidget())

        FigureCanvas.setSizePolicy(canvas,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(canvas)

        layout.addWidget(toolbar)
        layout.addWidget(canvas)

        return

    def draw(self):
        self.canvas.draw()
        return

class Spectrum():
    def __init__(self, x, real, filename=None, original_format=None):
        self.x = np.array(x)
        self.real = np.array(real)
        self.filename = filename
        self.original_format = original_format
        return

    @classmethod
    def open(cls, file, format = None):
        if format == None:
            if file.endswith('.asc'):
                format = 'Bruker TopSpin'
            else:
                format = 'SpinWorks'

        if format == 'Bruker TopSpin':
            datapoints = np.genfromtxt(file, skip_header = 2)
        elif format == 'SpinWorks':
            datapoints = np.genfromtxt(file, skip_header = 1)

        if datapoints[0][0] > datapoints[-1][0]:
            datapoints = np.flipud(datapoints)

        return Spectrum(datapoints[:, 0],
                        datapoints[:, 1],
                        filename=os.path.basename(file),
                        original_format=format)

    @classmethod
    def average_spectrum(cls, spectrums):
        return np.average(spectrums, axis = 0)

    @classmethod
    def similarity(cls, sp1, sp2):
        if len(sp1.x) == len(sp2.x):
            return Spectrum.similarity_samelen(sp1, sp2)
        elif len(sp1.x) > len(sp2.x):
            return Spectrum.similarity_difflen(sp2, sp1)
        else:
            return Spectrum.similarity_difflen(sp1, sp2)

    @classmethod
    def similarity_difflen(cls, sp_short, sp_long):
        sum1 = 0
        sum2 = 0
        intersection = 0
        for (x, real_short) in zip(sp_short.x, sp_short.real):
            real_long = sp_long.real_at(x)
            sum1 += real_short
            sum2 += real_long
            intersection += min(real_short, real_long)

        union = sum1 + sum2 - intersection
        #index = intersection / union
        index = intersection / 2

        return index

    @classmethod
    def similarity_samelen(cls, sp1, sp2):
        total = sp1.sum() + sp2.sum()
        diff_array = abs(sp1.real - sp2.real)
        diff = sum(diff_array)
        intersection = (total - diff)/2

        return intersection

        # union = total - intersection
        # index = intersection / union
        # return index

    @classmethod
    def average_spectrum(cls, sp_list):
        return Spectrum(np.copy(sp_list[0].x),
                        np.average(list(map(lambda sp: sp.real, sp_list)), axis = 0))

    def copy(self):
        return Spectrum(np.copy(self.x),
                        np.copy(self.real))

    def plot(self):
        plt.plot(self.x, self.real)
        return plt.show()

    @classmethod
    def plot_comparison(cls, sp1, sp2):
        difference = (sp1.real - sp2.real) * 8
        plt.plot(sp1.x, difference, label='difference*8')
        plt.plot(sp1.x, sp1.real + 0.0005, label=sp1.filename)
        plt.plot(sp1.x, -sp2.real - 0.0005, label=sp2.filename)
        plt.legend()
        return plt.show()

    def max(self):
        return np.max(self.real)

    def sum(self):
        return np.sum(self.real)

    def real_scaled(self):
        return self.real/self.max()

    def real_at(self, x):
        right = np.searchsorted(self.x, x)
        right_x = self.x[right]
        left_x = self.x[right - 1]
        right_val = self.real[right]
        left_val = self.real[right - 1]

        return left_val + (right_val - left_val) * (x - left_x) / (right_x - left_x)

    def remove_range(self, left, right):
        left_i = np.searchsorted(self.x, left)
        right_i = np.searchsorted(self.x, right, side='right')

        left_segment = self.real[:left_i]
        right_segment = self.real[right_i:]

        zeros = np.zeros(right_i - left_i)

        self.real = np.concatenate((left_segment, zeros, right_segment))

        return self

    def scale_max(self):
        self.real /= self.max()
        return self

    def scale_sum(self):
        self.real /= self.sum()
        return self

    def offset(self, n):
        zeros = np.zeros(abs(n))
        # 0 must match first clause since [0:] will get the whole array, [:0] won't
        if n >= 0:
            spectrum = Spectrum(self.x, np.concatenate((self.real[n:], zeros)))
        else:
            spectrum = Spectrum(self.x, np.concatenate((zeros, self.real[:n])))

        return spectrum

    def snip_edges(self, left = -0.5, right = 10):
        return self.remove_range(self.x[0], left).remove_range(right, self.x[-1])

    #def baseline(self, left = 180, right = 200):
    def baseline(self, left=10, right=12):
        left_i = np.searchsorted(self.x, left)
        right_i = np.searchsorted(self.x, right)
        baseline = np.average(self.real[left_i:right_i])

        self.real -= baseline

        return self

class SimilarityMatrix:
    def __init__(self, files, max_offset=10, transform=lambda sp: sp.snip_edges(left=-0.3, right=10.1).scale_sum()):
        self.max_offset = max_offset

        print('open')
        self.sp_list = list(map(lambda file: Spectrum.open(file), files))
        print('transform')
        self.sp_list = list(map(transform, self.sp_list))
        self.alignment_matrix = np.zeros((len(self.sp_list), 2 * max_offset + 1))
        self.best_alignments = np.zeros(len(self.sp_list))

        self.average_spectrum = self.find_average_spectrum()
        return

    @classmethod
    def open_standard_files(cls, max_offset=10):
        return cls(files, max_offset)

    def find_average_spectrum(self):
        return Spectrum.average_spectrum(self.sp_list)

    def find_best_alignments(self):
        for sp_i, sp in enumerate(self.sp_list):
            for off_i, off in enumerate(range(-self.max_offset, self.max_offset + 1)):
                self.alignment_matrix[sp_i, off_i] = Spectrum.similarity(sp.offset(off), self.average_spectrum)

        self.best_alignments = np.argmax(self.alignment_matrix, axis = 1) - self.max_offset
        return self.best_alignments

    def offset_aligned_spectrums(self):
        aligned_sp_list = []
        for sp, off in zip(self.sp_list, self.best_alignments):
            aligned_sp_list.append(sp.offset(off))
        return aligned_sp_list

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='*')

    parser.add_argument('--no-baseline-adjust', dest = 'baseline_adjust',
                        action = 'store_const', const = False, default = True)
    parser.add_argument('--no-edge-snip', dest = 'edge_snip',
                        action = 'store_const', const = False, default = True)
    parser.add_argument('--console', action = 'store_const', const = True, default = False)

    return parser.parse_args()

def console_execution(files = None):
    return

files = ["data/D2 1Hdump","data/DDOM1","data/DDOM2 1Hdump","data/DDOM3 1Hdump","data/DDOM4 1Hdump","data/JSAL1dump","data/js-alm1","data/JSALQ1","data/JSALQ2","data/JSALQ3","data/js-c1","data/JSC1Xdump","data/js-c2","data/JSC2Xdump","data/js-ca1","data/js-ca2","data/js-ca3","data/JSCA4dump","data/JSCA4X dump","data/JSCA5dump","data/JSCA5X dump","data/JSCA6dump","data/JSCHIA","data/js-d1","data/js-e1","data/js-g1","data/js-g2","data/js-g3","data/JSG4dump","data/js-g5","data/js-l1","data/js-l2","data/JSL3dump","data/js-l4","data/js-m1","data/js-m2","data/js-m3","data/js-mac1","data/js-ol1","data/JSOL2dump","data/js-ol3c","data/js-ol4","data/js-ol5c","data/js-s1","data/js-s2","data/JSS3dump","data/JSS4dump","data/js-s5","data/js-s6","data/js-wal1","data/OL6 1Hdump","data/OL7 1Hdump","data/SESAC 1Hdump","data/WAL21 1Hdump","data/WAL22 1Hdump"]

if __name__ == '__main__':
    parsed = parse_args(sys.argv)

    if parsed.console:
        #ret = console_execution(files = parsed.files, baseline_adjust = parsed.baseline_adjust, edge_snip = parsed.edge_snip)
        ret = console_execution(files = parsed.files)
        sys.exit(ret)
    else:
        print(sys.argv)
        print(parsed)
        app = QApplication(sys.argv)
        ex = App(files = parsed.files)
        #ex = App(files = parsed.files, baseline_adjust = parsed.baseline_adust, edge_snip = parsed.edge_snip)
        sys.exit(app.exec_())
