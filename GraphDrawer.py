from GoldbergWeighted import WFind_Densest_Subgraph, WFind_Density

import numpy as np
import matplotlib.pyplot as plt

def main():
    filepath = "rand_breakpoints.txt"
    file = open(filepath, "r")

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # Data for a three-dimensional line
    zline = np.linspace(0, 15, 1000)
    xline = np.sin(zline)
    yline = np.cos(zline)
    ax.plot3D(xline, yline, zline, 'gray')

    # Data for three-dimensional scattered points
    zdata = 15 * np.random.random(100)
    xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
    ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
    ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');
    plt.show()


    number_of_nodes = 499
    number_of_edges = 105
    answer_pred = WFind_Densest_Subgraph(number_of_nodes, number_of_edges, filepath)
    dens_pred = WFind_Density(answer_pred, filepath)
    print(dens_pred)
    print(answer_pred)
    file.close()

if __name__ == '__main__':
    main()