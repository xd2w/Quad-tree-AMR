import os
from subprocess import Popen

global text

text = "p "


def show(t=10):
    global text
    with open("temp.plt", "w") as f:
        filename = "out"

        preamble = f"""
        set xrange [0:1]
        #set yrange [-0.5:0.5]
        set yrange [0:1]

        set size ratio -1

        """
        f.write(preamble + text + f"\npause {t}\n")

    with Popen("gnuplot temp.plt", shell=True) as p:
        p.wait()


def save(filename):
    global text
    with open("temp.plt", "w") as f:
        filename = "out"

        preamble = f"""
        set xrange [0:1]
        #set yrange [-0.5:0.5]
        set yrange [0:1]

        set size ratio -1

        set term eps
        set output "{filename}.eps"

        """
        f.write(preamble + text + "\n")

    with Popen("gnuplot temp.plt", shell=True) as p:
        p.wait()


def plot_200(dir):
    global text
    text += f"'../{dir}/DATA/intf.200' w l, "


def plot_400(dir):
    global text
    text += f"'../{dir}/DATA/intf.400' w l, "


def plot_th200(dir):
    global text
    text += f"'../{dir}/DATA/thintf.200' w l, "


def plot_th400(dir):
    global text
    text += f"'../{dir}/DATA/thintf.000' w l, "


plot_200("Lv9_3")
plot_200("Lv10_3")
plot_200("Lv11_3")
plot_200("Lv12_3")
plot_th200("Lv12_3")
show()
