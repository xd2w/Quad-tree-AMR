#!/usr/bin/env python

import argparse

import numpy as np
from matplotlib import pyplot as plt

def main():
    par = argparse.ArgumentParser()
    par.add_argument('path')
    par.add_argument('-c', '--color')
    arg = par.parse_args()

    if arg.color is None:
        arg.color = 'r'

    with open(arg.path) as f:
        data = f.read()

    data = data.strip().splitlines()
    data = [list(map(float, datum.strip().split())) for datum in data]
    #data = [row if len(row) < 2 else continue for row in data]

   # print(data[0])
    
    x=[]
    y=[]

    for row in data:
        if len(row) < 2:
            plt.plot(x, y, arg.color)
            x = []
            y = []
            continue
        x.append(row[0])
        y.append(row[1])

    #x = np.array(x)
    #y = np.array(y)

    plt.plot(x, y, arg.color)
    plt.gca().set_aspect('equal')
    plt.show()

if __name__ == "__main__":
    main()
