import os
from subprocess import Popen

global text

text = "p "


def show(t=1, limits=False, equal_axes=True):
    global text
    with open("temp.plt", "w") as f:
        preamble = ""
        if equal_axes:
            preamble += """
            set size ratio -1
            """
        if limits:
            preamble += """
            set xrange [0:1]
            set yrange [0:1]
            """

        f.write(preamble + text + f"\npause {t}\n")

    text = "p "

    with Popen("gnuplot temp.plt", shell=True) as p:
        p.wait()


def save(filename="out"):
    global text
    with open("temp.plt", "w") as f:
        preamble = f"""
        set xrange [0:1]
        #set yrange [-0.5:0.5]
        set yrange [0:1]

        set size ratio -1

        set term eps
        set output "{filename}.eps"

        """
        f.write(preamble + text + "\n")

    text = "p "

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


def plot_intf(dir, n):
    global text
    text += f"'../{dir}/DATA/intf.{n:03d}' w l, "


def plot_mesh(dir, n):
    global text
    text += f"'../{dir}/DATA/mesh.{n:03d}' w l, "


def plot_kappa(dir, n):
    global text
    text += f"'../{dir}/DATA/kappa.{n:03d}' u 1:2, '../{dir}/DATA/kappa.{n:03d}' u 1:3 "


def plot_n_with_kappa(dir, n):
    global text
    text += f"'../{dir}/DATA/fgrd.{n:03d}' u 1:2:($3/($5+1e-50)*1e-3):($4/($5+1e-50)*1e-3) w vector"


# plot_400("Lv9_3")
# plot_400("Lv10_3")
# plot_400("Lv11_3")
# plot_400("Lv12_3")
# plot_th400("Lv12_3")
# show()

plot_n_with_kappa("HF", 0)
# plot_kappa("HF", 0)
show(limits=False, equal_axes=False)

plot_n_with_kappa("CSF", 0)
# plot_kappa("HF", 0)
show(limits=False, equal_axes=False)

# plot_mesh("Lv9_3", 20)
# plot_intf("Lv9_3", 20)
# # show()
# save("kappa_example")

# plot_mesh("Lv9_Uniform", 20)
# plot_intf("Lv9_Uniform", 20)
# # show()
# save("uniform_example")
