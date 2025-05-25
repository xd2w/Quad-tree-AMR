import glob
import re

# from collections import Iterator
from collections.abc import Iterable

import numpy as np
from matplotlib import pyplot as plt

dirs = glob.glob("../**/DATA", recursive=True)


def total_computation(dir):
    data = np.loadtxt(f"../{dir}/DATA/runtime.data", delimiter=",", skiprows=4)
    # print(data.shape)
    iterations = data[:, 0]
    runtimes = data[:, 1]
    dts = data[:, 2]
    cells = data[:, 3]

    # np.trapezoid(iterations, cells, initial=0)
    print(f"{dir} \t: {np.sum(cells):.3g}")


def cell_evolution(dir):
    data = np.loadtxt(f"../{dir}/DATA/runtime.data", delimiter=",", skiprows=4)
    # print(data.shape)
    iteration = data[:, 0]
    runtime = data[:, 1]
    dts = data[:, 2]
    cells = data[:, 3]

    # plt.plot(runtime, cells * iteration[-1] / 4, label=dir)
    plt.plot(runtime, cells / dts, label=dir)

    # np.trapezoid(iterations, cells, initial=0)
    print(f"{dir} \t: {np.sum(cells):.3g}")


def area_error(dir, fig_name=None):
    data = np.loadtxt(f"../{dir}/DATA/intf.400")
    r = 0.2
    xc = 0.5
    yc = 0.75
    origin = np.array([xc, yc])
    # print(data)
    area = 0.0
    for i in range(0, data.shape[0], 2):
        r1 = data[i]
        r2 = data[i + 1]

        theta1 = np.arctan2(r1[1] - yc, r1[0] - xc)
        theta2 = np.arctan2(r2[1] - yc, r2[0] - xc)

        l1 = np.array([r * np.cos(theta1), r * np.sin(theta1)]) + origin
        l2 = np.array([r * np.cos(theta2), r * np.sin(theta2)]) + origin
        plt.fill(
            [l1[0], l2[0], r2[0], r1[0]], [l1[1], l2[1], r2[1], r1[1]], "r", alpha=1
        )
        plt.plot([r1[0], r2[0]], [r1[1], r2[1]], "b-")
        plt.plot([l1[0], l2[0]], [l1[1], l2[1]], "k-")

        area += 0.5 * np.linalg.norm(np.cross(l2 - l1, l2 - r2)) + 0.5 * np.linalg.norm(
            np.cross(r2 - r1, l1 - r1)
        )
    # area = np.pi * (r2**2 - r1**2)
    print(f"{dir} : Error Area: {area/(np.pi * r**2):.4g}%")
    plt.axis("equal")
    if fig_name is not None:
        plt.savefig(f"Figs/{fig_name}.eps")
        plt.cla()
    else:
        plt.show()


def plot_400(dir, label=None, color=None):

    if label is None:
        label = dirs

    if color is None:
        color = "r"
    elif isinstance(color, Iterable):
        color = next(color)

    data = np.loadtxt(f"../{dir}/DATA/intf.400")

    # r = 0.2
    # xc = 0.5
    # yc = 0.75
    # origin = np.array([xc, yc])
    for i in range(0, data.shape[0], 2):
        r1 = data[i]
        r2 = data[i + 1]

        plt.plot(
            [r1[0], r2[0]], [r1[1], r2[1]], color=color, label=label if i == 0 else None
        )


def plot_200(dir, label=None, color=None):

    if label is None:
        label = dirs

    if color is None:
        color = "r"
    elif isinstance(color, Iterable):
        color = next(color)

    data = np.loadtxt(f"../{dir}/DATA/intf.200")

    # r = 0.2
    # xc = 0.5
    # yc = 0.75
    # origin = np.array([xc, yc])
    for i in range(0, data.shape[0], 2):
        r1 = data[i]
        r2 = data[i + 1]

        plt.plot(
            [r1[0], r2[0]], [r1[1], r2[1]], color=color, label=label if i == 0 else None
        )


def plot_th200(dir, label=None, color=None):

    if label is None:
        label = dirs

    if color is None:
        color = "r"
    elif isinstance(color, str):
        pass
    elif isinstance(color, Iterable):
        color = next(color)

    data = np.loadtxt(f"../{dir}/DATA/thintf.200")
    plt.plot(data[:, 0], data[:, 1], color, label=label)


def plot_circle(xc, yc, r, **kwargs):
    theta = np.linspace(0, 2 * np.pi, 100)
    x = xc + r * np.cos(theta)
    y = yc + r * np.sin(theta)
    plt.plot(x, y, **kwargs)


def str2time(time_str):

    pattern = r"(?:(\d+)m)?([\d.]+)s"
    match = re.match(pattern, time_str)

    if not match:
        raise ValueError(f"Invalid time format: {time_str}")

    minutes = int(match.group(1)) if match.group(1) else 0
    seconds = float(match.group(2))

    return minutes * 60 + seconds


def cpu_time(dir):
    print(f"{dir}:")
    with open(f"../{dir}/time.txt", "r") as f:
        data = f.read()
    data = data.strip()
    data = data.splitlines()
    data = [datum.split() for datum in data]
    data = {key: str2time(val) for key, val in data}
    print(f"{data["sys"] + data["user"]:.3g}")
    print()


def plot_step_drop(n=9):
    color = iter(["r", "g", "b", "c", "m", "y"])
    plot_400(f"Lv{n-1}_Uniform", color=color, label=f"Uniform {n-1}")
    plot_400(f"Lv{n}_3", color=color, label=f"Kappa {n}")
    plot_400(f"Lv{n}_Uniform", color=color, label=f"Uniform {n}")
    # plot_400("Lv12_3", color=color, label="Kappa 12")
    plot_circle(0.5, 0.75, 0.2, color="k", linestyle="--", label="Target Circle")
    plt.legend()
    plt.show()


def plot_all_400_uni():
    color = iter(["r", "g", "b", "c", "m", "y"])
    plot_400("Lv7_Uniform", color=color, label="Uniform 7")
    plot_400("Lv8_Uniform", color=color, label="Uniform 8")
    plot_400("Lv9_Uniform", color=color, label="Uniform 9")
    plot_400("Lv10_Uniform", color=color, label="Uniform 10")
    plot_circle(0.5, 0.75, 0.2, color="k", linestyle="--", label="Target Circle")
    plt.legend()
    plt.show()


def plot_all_400_kappa():
    color = iter(["r", "g", "b", "c", "m", "y"])
    plot_400("Lv9_3", color=color, label="Kappa 9")
    plot_400("Lv10_3", color=color, label="Kappa 10")
    plot_400("Lv11_3", color=color, label="Kappa 11")
    # plot_400("Lv12_3", color=color, label="Kappa 12")
    plot_circle(0.5, 0.75, 0.2, color="k", linestyle="--", label="Target Circle")
    plt.legend()
    plt.show()


def plot_all_200_uni():
    color = iter(["r", "g", "b", "c", "m", "y"])
    plot_th200("Lv12_3", color="k", label="Target")
    plot_200("Lv7_Uniform", color=color, label="Uniform 7")
    plot_200("Lv8_Uniform", color=color, label="Uniform 8")
    plot_200("Lv9_Uniform", color=color, label="Uniform 9")
    plot_200("Lv10_Uniform", color=color, label="Uniform 10")
    plt.legend(loc="upper right")
    plt.show()


def plot_all_200_kappa():
    color = iter(["r", "g", "b", "c", "m", "y"])
    plot_th200("Lv12_3", color="k", label="Target")
    plot_200("Lv9_3", color=color, label="Kappa 9")
    plot_200("Lv10_3", color=color, label="Kappa 10")
    plot_200("Lv11_3", color=color, label="Kappa 11")
    plot_200("Lv12_3", color=color, label="Kappa 12")
    plt.legend(loc="upper right")
    plt.show()


if __name__ == "__main__":
    # total_computation("Lv9_3")
    # total_computation("Lv8_Uniform")

    # area_error("Lv7_Uniform", "uni_7")
    # area_error("Lv8_Uniform", "uni_8")
    # area_error("Lv9_Uniform", "uni_9")
    # area_error("Lv10_Uniform", "uni_10")

    # area_error("Lv9_3", "kappa_9")
    # area_error("Lv10_3", "kappa_10")
    # area_error("Lv11_3", "kappa_11")
    # area_error("Lv12_3", "kappa_12")

    # cpu_time("Lv7_Uniform")
    # cpu_time("Lv8_Uniform")
    # cpu_time("Lv9_Uniform")
    # cpu_time("Lv10_Uniform")

    # cpu_time("Lv9_3")
    # cpu_time("Lv10_3")
    # cpu_time("Lv11_3")
    # cpu_time("Lv12_3")

    # total_computation("Lv7_Uniform")
    # total_computation("Lv8_Uniform")
    # total_computation("Lv9_Uniform")
    # total_computation("Lv10_Uniform")

    # total_computation("Lv9_3")
    # total_computation("Lv10_3")
    # total_computation("Lv11_3")
    # total_computation("Lv12_3")

    # plot_all_400_uni()
    # plot_all_400_kappa()

    # plot_step_drop(9)
    # plot_step_drop(10)

    # plot_all_200_uni()
    # plot_all_200_kappa()

    cell_evolution("Lv7_Uniform")
    cell_evolution("Lv8_Uniform")
    cell_evolution("Lv9_Uniform")
    cell_evolution("Lv10_Uniform")

    cell_evolution("Lv9_3")
    cell_evolution("Lv10_3")
    cell_evolution("Lv11_3")
    # cell_evolution("Lv12_3")

    plt.legend()
    plt.show()
