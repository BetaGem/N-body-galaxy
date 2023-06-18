import numpy as np

# spherical symmetric profiles

def const(r, args):
    '''
    args[0]: a constant
    '''
    return args[0]


def Gaussian(r, args):
    '''
    args[0]: normalization
    args[1]: std
    '''
    return args[0] / args[1] / np.exp(r ** 2 / 2 / args[1] ** 2)


def power(r, args):
    '''
    args[0]: normalization
    args[1]: characteristic radius
    args[2]: power law index
    '''
    return args[0] * (args[1] / r) ** args[2]


def double_power(r, args):
    '''
    args[0]: normalization
    args[1]: characteristic radius
    args[2]: inner power law index
    args[3]: outer power law index
    '''
    return args[0] / (r / args[1]) ** args[2] / (1 + r / args[1]) ** (args[3] - args[2])