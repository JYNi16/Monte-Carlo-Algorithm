import matplotlib.pylab as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np
from xymodel_mc import *

def visualise(sim):
    X = np.arange(sim.A.size).reshape(sim.A.shape) % sim.A.shape[0]
    Y = (np.arange(sim.A.size).reshape(sim.A.shape) % sim.A.shape[1]).T
    #X, Y = np.meshgrid((0, sim.A[0], .2), np.arange(0, sim.A[1], .2))

    U = cos(sim.A)
    V = sin(sim.A)

    fig, ax = plt.subplots()

    # rects = []
    #
    # # create rectangles for vortex/antivortex determination
    # for i in range(sim.V.shape[0]):
    #     for j in range(sim.V.shape[1]):
    #         rect = patches.Rectangle(xy=(i, j), height=1, width=1)
    #         rects.append(rect)
    # rects = PatchCollection(rects)
    #
    # # Set colors for the rectangles
    # col = 'RdBu'
    # r_cmap = plt.get_cmap(col)
    # r_cmap_r = plt.get_cmap(col + "_r")  # eto kostil' =)
    # rects.set_cmap(r_cmap)
    # rects.set_clim(vmin=-1, vmax=1)
    #
    # rects.set_animated(True)
    # rects.set_array(sim.V.flatten('F') / 2)
    # ax.add_collection(rects)
    #
    # # create legend
    # legend_boxes = [patches.Patch(facecolor=r_cmap(0.7), label='Antiortex'),
    #                 patches.Patch(facecolor=r_cmap_r(0.7), label='Vortex')]
    # ax.legend(handles=legend_boxes)

    # build an initial quiver plot
    #q = ax.quiver(X, Y, U, V, pivot='mid', cmap=plt.cm.get_cmap('hsv'), units='inches', scale=4)
    q = ax.quiver(X, Y, U, V, units='width')
    qk = ax.quiverkey(q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',coordinates='figure')
    #fig.colorbar(q, label='Angles (2 pi)')
    # ax.set_xlim(-1, sim.A.shape[0])
    # ax.set_ylim(-1, sim.A.shape[1])
    #qk.set_UVC(U, V, C=sim.A)

    plt.show()

    return q, fig


if __name__=="__main__":
    sim = XYMetropolis((50,50),
                 beta=1,
                 J=10,
                 random_state=5,
                 initial_state='hot')
    sim.simulate(1000000)
    visualise(sim)
