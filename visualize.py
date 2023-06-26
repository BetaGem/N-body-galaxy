import numpy as np
import matplotlib.pyplot as plt
from photutils.segmentation import make_2dgaussian_kernel
from astropy.convolution import convolve

from mass_assign import CIC_3D
from params import *

def plot_all(part, dens, Phi, acc, save=False):
    '''
    Description: plot the result.
    '''
    plt.clf()
    plt.figure(1, figsize=(16,4))

    plt.subplot(141, title="particles")
    plt.scatter(part[:,ix], part[:,iy], c="k", s=.1, alpha=.3)
    plt.xlim(0, Ng); plt.ylim(0, Ng)
    plt.subplot(142, title="density field")
    plt.imshow(np.arcsinh(np.sum(dens, axis=0)), origin="lower")
    plt.subplot(143, title="gravitational potential")
    plt.imshow(np.sum(Phi, axis=0), origin="lower", cmap="jet")
    plt.clim(np.min(np.sum(Phi, axis=0)), np.max(np.sum(Phi, axis=0)))
    plt.subplot(144, title="acceleration field")
    plt.imshow(acc[1][z] ** 2 + acc[2][z] ** 2, origin="lower", cmap="jet")
    plt.clim(np.min(acc[1] ** 2 + acc[2] ** 2) ,np.max(acc[1] ** 2 + acc[2] ** 2))
    plt.pause(0.01)
    if save: plt.save(save)
    
def save_data(data, path="."):
    '''
    Description: save particle info to files
    '''
    np.save(path, data)


def plot_file(path, slice_num):
    '''
    Description: generate animations from files
    '''
    plt.figure(1, figsize=(6,6))
    
    for i in range(slice_num):
        # load particles
        part = np.load(path + str(i) + ".npy")
        # animation
        plt.clf()
        # plt.subplot(121)
        # plt.xlim(0, Ng); plt.ylim(0, Ng)
        # plt.scatter(part[:,ix], part[:,iy], c="k", s=.1, alpha=.2)
        #
        # plt.subplot(122)
        dens = CIC_3D(part)
        kernel = make_2dgaussian_kernel(fwhm=5, size=15)
        intens = convolve(np.sum(dens, axis=0) + np.random.randn(512, 512)*.01, kernel)
        intens[intens < 0] = 0
        plt.imshow(np.log10(intens + .02), cmap="gray", origin="lower")
        plt.pause(0.01)

if __name__ == "__main__":

    plot_file("./data/", slice_num=int(t_tot/dt))
