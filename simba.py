#!/usr/bin/env python3

import os
import sys
import random
import argparse
import re
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from typing import List


class Spectrum:
    def __init__(self, x, real, filename=None, original_format=None):
        self.x = np.array(x)
        self.real = np.array(real)
        self.filename = filename
        self.original_format = original_format
        return

    @classmethod
    def open(cls, file, format=None):
        if format is None:
            if file.endswith(".asc"):
                format = "Bruker TopSpin"
            else:
                format = "SpinWorks"

        if format == "Bruker TopSpin":
            datapoints = np.genfromtxt(file, skip_header=1)
            datapoints = np.flipud(datapoints)
            x = datapoints[:, 0]
            real = datapoints[:, 1]
        elif format == "SpinWorks":
            with open(file, 'r', encoding='latin1', newline='\r\n') as f:
                header = ""
                for _ in range(2):
                    header += f.readline()

                datapoints = np.flipud(np.genfromtxt(f))

                x = datapoints[:, 0]
                real = datapoints[:, 1]
        elif format == "TopSpin2":
            with open(file, 'r', encoding='latin1') as f:
                header = ""
                for _ in range(10):
                    header += f.readline()

                left = float(re.search("LEFT = (-?[0-9.]+)", header).group(1))
                right = float(re.search("RIGHT = (-?[0-9.]+)", header).group(1))
                size = int(re.search("SIZE = ([0-9]+)", header).group(1))
                x = np.linspace(right, left, size)
                real = np.flipud(np.genfromtxt(f))
        elif format == "TSV":
            datapoints = np.genfromtxt(file, skip_header=0)
            datapoints = np.flipud(datapoints)
            x = datapoints[:, 0]
            real = datapoints[:, 1]
        elif format == "ACD":
            with open(file, 'r', encoding='latin1') as f:
                header = ""

                line = f.readline()
                while line.strip() != "":
                    header += line
                    line = f.readline()

                datapoints = np.flipud(np.genfromtxt(f))
                x = datapoints[:, 0]
                real = datapoints[:, 1]

        return Spectrum(
            x,
            real,
            filename=os.path.basename(file),
            original_format=format,
        )

    @classmethod
    def open_spectra(cls, files, format=None):
        return list(map(lambda file: Spectrum.open(file, format), files))

    # @classmethod
    # def average_spectrum(cls, spectra):
    #     real_axes = list(map(lambda sp: sp.real, spectra))
    #     return Spectrum(spectra[0].x, np.average(real_axes, axis=0))

    @classmethod
    def average_spectrum(cls, sp_list):
        return Spectrum(
            np.copy(sp_list[0].x),
            np.average(list(map(lambda sp: sp.real, sp_list)), axis=0),
        )

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
        # index = intersection / union
        index = intersection / 2

        return index

    @classmethod
    def similarity_samelen(cls, sp1, sp2):
        total = sp1.sum() + sp2.sum()
        diff_array = abs(sp1.real - sp2.real)
        diff = sum(diff_array)
        intersection = (total - diff) / 2

        return intersection

        # union = total - intersection
        # index = intersection / union
        # return index

    def copy(self):
        return Spectrum(np.copy(self.x), np.copy(self.real))

    def plot(self):
        plt.plot(self.x, self.real)
        return plt.show()

    @classmethod
    def plot_comparison(cls, sp1, sp2):
        difference = (sp1.real - sp2.real) * 8
        plt.plot(sp1.x, difference, label="difference*8")
        plt.plot(sp1.x, sp1.real + 0.0005, label=sp1.filename)
        plt.plot(sp1.x, -sp2.real - 0.0005, label=sp2.filename)
        plt.legend()
        return plt.show()

    def max(self):
        return np.max(self.real)

    def sum(self):
        return np.sum(self.real)

    def real_scaled(self):
        return self.real / self.max()

    def real_at(self, x):
        right = np.searchsorted(self.x, x)
        right_x = self.x[right]
        left_x = self.x[right - 1]
        right_val = self.real[right]
        left_val = self.real[right - 1]

        return left_val + (right_val - left_val) * (x - left_x) / (right_x - left_x)

    def remove_range(self, left, right):
        left_i = np.searchsorted(self.x, left)
        right_i = np.searchsorted(self.x, right, side="right")

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
            spectrum = Spectrum(
                self.x, np.concatenate((self.real[n:], zeros)), filename=self.filename
            )
        else:
            spectrum = Spectrum(
                self.x, np.concatenate((zeros, self.real[:n])), filename=self.filename
            )

        return spectrum

    def snip_edges(self, left=-0.5, right=10):
        return self.remove_range(self.x[0], left).remove_range(right, self.x[-1])

    # def baseline(self, left = 180, right = 200):
    def baseline(self, left=10, right=12):
        left_i = np.searchsorted(self.x, left)
        right_i = np.searchsorted(self.x, right)
        baseline = np.average(self.real[left_i:right_i])

        self.real -= baseline

        return self

    def rolling_average(self, n):
        conv_filter = np.repeat(1 / n, n)
        ltail = n // 2
        rtail = (n - 1) // 2
        self.real = np.convolve(self.real, conv_filter)[ltail:-rtail]

        return self

    def cancel_negative(self):
        sorted_real = np.sort(self.real)
        cumsum = np.cumsum(sorted_real)

        threshold_index = np.searchsorted(cumsum, 0)
        threshold = sorted_real[threshold_index - 1]

        print(threshold)

        self.real[self.real < threshold] = 0

        return self

    def output(self, format="SpinWorks"):

        return


class SimilarityMatrix:
    def __init__(self, sp_list, max_offset=10):
        self.max_offset = max_offset

        self.sp_list = sp_list
        self.alignment_matrix = np.zeros((len(self.sp_list), 2 * max_offset + 1))
        self.similarity_matrix = np.zeros((len(self.sp_list), len(self.sp_list)))
        self.best_alignments = np.zeros(len(self.sp_list))

        self.average_spectrum = self.find_average_spectrum()
        return

    def find_average_spectrum(self):
        return Spectrum.average_spectrum(self.sp_list)

    def find_best_alignments(self):
        for sp_i, sp in enumerate(self.sp_list):
            for off_i, off in enumerate(range(-self.max_offset, self.max_offset + 1)):
                self.alignment_matrix[sp_i, off_i] = Spectrum.similarity(
                    sp.offset(off), self.average_spectrum
                )

        self.best_alignments = (
            np.argmax(self.alignment_matrix, axis=1) - self.max_offset
        )
        return self.best_alignments

    def offset_aligned_spectrums(self):
        aligned_sp_list = []
        for sp, off in zip(self.sp_list, self.best_alignments):
            aligned_sp_list.append(sp.offset(off))
        return aligned_sp_list

    def populate_matrix(self):
        for i in range(0, len(self.sp_list)):
            for j in range(0, len(self.sp_list)):
                if i < j:
                    self.similarity_matrix[i, j] = Spectrum.similarity(
                        self.sp_list[i], self.sp_list[j]
                    )
                elif i == j:
                    self.similarity_matrix[i, j] = 1
                else:
                    self.similarity_matrix[i, j] = self.similarity_matrix[j, i]
                    # self.similarity_matrix[i, j] = self.similarity_matrix[j, i]
        return self.similarity_matrix


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*")

    parser.add_argument(
        "--no-baseline-adjust",
        dest="baseline_adjust",
        action="store_const",
        const=False,
        default=True,
    )
    parser.add_argument(
        "--no-edge-snip",
        dest="edge_snip",
        action="store_const",
        const=False,
        default=True,
    )
    parser.add_argument("--console", action="store_const", const=True, default=False)

    return parser.parse_args()


from dataclasses import dataclass
from os import listdir
from os.path import isfile, join


@dataclass
class PipelineStep:
    name: str
    options: dict = None

    def __repr__(self):
        return f"Step({self.name}, {self.options})"


class Arguments:
    def __init__(self, args):
        self.files = []
        self.pipeline = []
        self.format = None

        args = iter(args)
        next(args)

        for arg in args:
            if arg[0] != "-":
                self.files.append(arg)

            elif arg == "--directory":
                directory = next(args)
                files = [
                    join(directory, f)
                    for f in listdir(directory)
                    if isfile(join(directory, f))
                ]
                self.files.extend(files)

            elif arg == "--baseline-adjust":
                limits = next(args).split(",", maxsplit=1)
                limits = list(map(lambda s: float(s), limits))
                self.pipeline.append(PipelineStep("baseline_adjust"))

            elif arg == "--snip-edges":
                limits = next(args).split(",", maxsplit=1)
                limits = list(map(lambda s: float(s), limits))
                self.pipeline.append(
                    PipelineStep("snip_edges", {"left": limits[0], "right": limits[1]})
                )

            elif arg == "--remove-range":
                limits = next(args).split(",", maxsplit=1)
                limits = list(map(lambda s: float(s), limits))
                self.pipeline.append(
                    PipelineStep(
                        "remove_range", {"left": limits[0], "right": limits[1]}
                    )
                )

            elif arg == "--scale-sum":
                self.pipeline.append(PipelineStep("scale_sum"))

            elif arg == "--rolling-average":
                n = int(next(args))
                self.pipeline.append(PipelineStep("rolling_average", {"n": n}))

            elif arg == "--cancel-negative":
                self.pipeline.append(PipelineStep("cancel_negative"))

            elif arg == "--similarity-align":
                offset = int(next(args))
                self.pipeline.append(
                    PipelineStep("similarity_align", {"max_offset": offset})
                )

            elif arg == "--similarity-matrix":
                self.pipeline.append(PipelineStep("similarity_matrix"))

            elif arg == "--format":
                self.format = next(args)

            else:
                print(f"Error: Unrecognized argument {arg}")
                sys.exit(1)

    def __repr__(self):
        return f"Files: {self.files}\nPipeline: {self.pipeline}"


class PipelineExecutor:
    def __init__(self, pipeline: List[PipelineStep], files: List[str], format: str):
        self.files = files
        self.format = format
        self.pipeline = pipeline

    def run(self):
        self.sp_list = Spectrum.open_spectra(self.files, format=self.format)

        for step in self.pipeline:
            self.run_step(step)
        return

    def run_step(self, step):
        print(f"{step.name}")
        if step.name == "baseline_adjust":
            for sp in self.sp_list:
                sp.baseline(step.options["left"], step.options["right"])

        elif step.name == "snip_edges":
            for sp in self.sp_list:
                sp.snip_edges(step.options["left"], step.options["right"])

        elif step.name == "remove_range":
            for sp in self.sp_list:
                sp.remove_range(step.options["left"], step.options["right"])

        elif step.name == "scale_sum":
            for sp in self.sp_list:
                sp.scale_sum()

        elif step.name == "rolling_average":
            for sp in self.sp_list:
                sp.rolling_average(step.options["n"])

        elif step.name == "cancel_negative":
            for sp in self.sp_list:
                sp.cancel_negative()

        elif step.name == "similarity_align":
            self.sm = SimilarityMatrix(self.sp_list, step.options["max_offset"])
            self.sm.find_best_alignments()
            self.old_sp_list = self.sp_list
            self.sm.sp_list = self.sm.offset_aligned_spectrums()
            self.sp_list = self.sm.sp_list

        elif step.name == "similarity_matrix":
            self.sm.populate_matrix()

        # elif step.name == 'out_matrix':
        #     filenames = map(lambda sp: sp.filename, self.sp_list)
        #     header = ','.join(filenames)
        #     np.savetxt(step.options['filename'], delimiter=',', header=header)

        else:
            print(f"Error: Unrecognized step {step.name}")
            sys.exit(1)


def output_matrix(filename, sm):
    filenames = map(lambda sp: sp.filename, sm.sp_list)
    header = ",".join(filenames)
    np.savetxt(
        filename, sm.similarity_matrix, delimiter=",", comments="", header=header
    )


def output_alignments(filename, sm):
    filenames = map(lambda sp: sp.filename, sm.sp_list)
    max_alignment = (len(sm.alignment_matrix[0]) - 1) / 2
    alignments = range(int(-max_alignment), int(max_alignment + 1))
    columns = ["best alignment"] + list(alignments)
    best_alignments = np.transpose([sm.best_alignments])
    table = np.concatenate([best_alignments, sm.alignment_matrix], axis=1)
    df = pd.DataFrame(table, index=filenames, columns=columns)
    df.to_csv(filename)


def output_aligned(sm, filename):
    reals = list(map(lambda sp: sp.real, sm.sp_list))
    filenames = list(map(lambda sp: sp.filename, sm.sp_list))
    df = pd.DataFrame(reals, index=filenames, columns=sm.sp_list[0].x)
    df.transpose().to_csv(filename)


def main(args):
    parsed = Arguments(sys.argv)
    print(parsed)

    executor = PipelineExecutor(parsed.pipeline, parsed.files, parsed.format)

    executor.run()

    output_matrix("out_matrix", executor.sm)
    output_alignments("out_align", executor.sm)
    output_aligned(executor.sm, "out_aligned")


if __name__ == '__main__':
    main(sys.argv)
