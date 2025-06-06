import os
import shutil
import glob
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


def save(filename="out", limits=False, equal_axes=True):
    global text
    with open("temp.plt", "w") as f:
        preamble = f"""
        set term eps
        set output "{filename}.eps"
        """
        if equal_axes:
            preamble += """
            set size ratio -1
            """
        if limits:
            preamble += """
            set xrange [0.7:0.9]
            set yrange [0.7:0.9]
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
    text += f"'../{dir}/DATA/intf.400' w l notitle, "


def plot_th200(dir):
    global text
    text += f"'../{dir}/DATA/thintf.200' w l notitle, "


def plot_th400(dir):
    global text
    text += f"'../{dir}/DATA/thintf.000' w l notitle, "


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
    text += f"'../{dir}/DATA/intf.{n:03d}' w l notitle, "
    text += f"'../{dir}/DATA/fgrd.{n:03d}' u 1:2:($3*($5)*4e-5):($4*($5)*4e-5) w vector notitle"


def plot_wl(path):
    global text
    text += f"'{path}' w l notitle, "


# plot_400("Lv9_3")
# plot_400("Lv10_3")
# plot_400("Lv11_3")
# plot_400("Lv12_3")
# plot_th400("Lv12_3")
# show()

# plot_mesh("../HF", 0)
# plot_n_with_kappa("../HF", 0)
# plot_kappa("HF", 0)
# show(limits=False, equal_axes=True)
# save("HF_ellipse_fgrd")

# plot_mesh("../CSF", 0)
# plot_n_with_kappa("../CSF", 0)
# plot_kappa("HF", 0)
# show(limits=False, equal_axes=True)
# save("CSF_ellipse_fgrid")

# plot_mesh("Lv9_3", 20)
# plot_intf("Lv9_3", 20)
# # show()
# save("kappa_example")

# plot_mesh("Lv9_Uniform", 20)
# plot_intf("Lv9_Uniform", 20)
# # show()
# save("uniform_example")

dirs = glob.glob("../**/DATA/*.eps")

# for path in dirs:
#     plot_th400("Lv9_3")
#     plot_wl(path)
#     save(path, limits=True)

for path in dirs:
    name = path.replace("/DATA/", "_")
    name = "all_figs/" + name[3:-4].replace(".", "_") + ".eps"
    shutil.copyfile(path, name)
    print(name)

# print(dirs)
