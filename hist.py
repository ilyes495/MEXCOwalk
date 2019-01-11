from utils import *
import numpy as np
import matplotlib.pyplot as plt

def plot_histogram_list(l, b):
    plt.hist(l, bins=b)
    plt.show()
