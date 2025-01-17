#!/usr/bin/env python

import argparse

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


def each_file():
    par = argparse.ArgumentParser()
    par.add_argument("path")
    par.add_argument("-c", "--color")
    arg = par.parse_args()

    if arg.color is None:
        arg.color = "r"

    with open(arg.path) as f:
        data = f.read()

    data = data.strip().splitlines()
    data = [list(map(float, datum.strip().split())) for datum in data]
    # data = [row if len(row) < 2 else continue for row in data]

    # print(data[0])

    x = []
    y = []

    for row in data:
        if len(row) < 2:
            plt.plot(x, y, arg.color)
            x = []
            y = []
            continue
        x.append(row[0])
        y.append(row[1])

    # x = np.array(x)
    # y = np.array(y)

    plt.plot(x, y, arg.color)
    plt.gca().set_aspect("equal")
    plt.show()


def show_file(i, file, lines=None, color="r"):
    try:
        with open(f"DATA/{file}.{(i//100)%10}{(i//10)%10}{i%10}") as f:
            data = f.read()
    except FileNotFoundError:
        print("no more files found")
        exit(0)

    data = data.strip().splitlines()
    data = [list(map(float, datum.strip().split())) for datum in data]
    # data = [row if len(row) < 2 else continue for row in data]

    # print(data[0])

    x = []
    y = []
    if lines is None:
        global glines
        lines = glines

    lines[0].remove()

    for row in data:
        if len(row) < 2:
            plt.plot(x, y, color)
            x = []
            y = []
            continue
        x.append(row[0])
        y.append(row[1])

        # x = np.array(x)
        # y = np.array(y)
    plt.gca().set_aspect("equal")
    lines = plt.plot(x, y, color)
    if lines is not None:
        return lines
    else:
        glines = lines


def show_multiple_files(i, files):
    for j in range(len(files)):
        multi_lines[j] = show_file(
            i, files[j], lines=multi_lines[j], color=["r", "b", "g", "y"][j]
        )


if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument("files", nargs="+", choices=["hilb", "intf", "sfcv", "mesh"])
    par.add_argument("-c", "--color")
    arg = par.parse_args()

    if arg.color is None:
        arg.color = "r"

    global glines
    global multi_lines

    glines = plt.plot([], [])
    multi_lines = [plt.plot([], []) for _ in arg.files]

    plt.gcf().set_size_inches(10, 10)
    plt.gca().set_aspect("equal")
    plt.gca().set_xlim(-1, 3)
    plt.gca().set_ylim(-1, 3)

    ani = FuncAnimation(plt.gcf(), show_multiple_files, fargs=(arg.files,), frames=1000)
    plt.show()
