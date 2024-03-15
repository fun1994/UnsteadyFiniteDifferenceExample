# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:13:22 2024

@author: HFKJ059
"""

import numpy as np
from matplotlib import pyplot as plt


def read_1d(filename):
    with open("./data/" + filename + ".txt", "r") as file:
        data = file.read()
    data = data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    data = np.array(data)
    return data

def read():
    x = read_1d("x")
    phi1 = read_1d("phi1")
    phi2 = read_1d("phi2")
    phi3 = read_1d("phi3")
    phi4 = read_1d("phi4")
    return x, phi1, phi2, phi3, phi4

def plot(x, phi1, phi2, phi3, phi4):
    plt.plot(x, phi1, label="explicit Euler")
    plt.plot(x, phi2, label="implicit Euler")
    plt.plot(x, phi3, label="Crank-Nicolson")
    plt.plot(x, phi4, label="three-time-level")
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("phi")
    plt.show()

def main():
    x, phi1, phi2, phi3, phi4 = read()
    plot(x, phi1, phi2, phi3, phi4)


main()
