import numpy as np
import matplotlib.pyplot as plt
from photutils.segmentation import make_2dgaussian_kernel
from astropy.convolution import convolve

from params import *

def plot_all(part, dens, Phi, acc):
    '''
    Description: plot the result.
    '''
    plt.clf()
    plt.figure(1, figsize=(8,4))

    #plt.subplot(141, title="particles")
    #plt.scatter(part[:,ix], part[:,iy], c="k", s=.1, alpha=.3)
    #plt.xlim(0, Ng); plt.ylim(0, Ng)
    plt.subplot(121, title="density field")
    plt.imshow(np.arcsinh(np.sum(dens, axis=0)), origin="lower")
    
    plt.subplot(122, title="gravitational potential")
    plt.imshow(np.sum(Phi, axis=0), origin="lower", cmap="jet")
    plt.clim(np.min(np.sum(Phi, axis=0)), np.max(np.sum(Phi, axis=0)))
    #plt.subplot(144, title="acceleration field")
    #plt.imshow(acc[1][z] ** 2 + acc[2][z] ** 2, origin="lower", cmap="jet")
    #plt.clim(np.min(acc[1] ** 2 + acc[2] ** 2) ,np.max(acc[1] ** 2 + acc[2] ** 2))
    plt.pause(0.01)