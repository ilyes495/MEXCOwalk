from utils import *
import numpy as np
import matplotlib.pyplot as plt

def plot_histogram_list(l, b):
    plt.hist(l, bins=b)
    plt.show()

l = []
with open("mutex_neighbourhood.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split()
        l.append(float(line[2]))

print(l)
plot_histogram_list(l, 100)
