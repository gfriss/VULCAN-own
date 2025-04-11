import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from cycler import cycler

def reset_plt(ticksize, fontsize, fxsize, fysize, grid = True):
    #plt.style.use('seaborn-white') # was seaborn-v0_8-white
    plt.rcParams['xtick.labelsize'] = ticksize
    plt.rcParams['ytick.labelsize'] = ticksize
    plt.rcParams['xtick.major.width'] = 1.25
    plt.rcParams['ytick.major.width'] = 1.25
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['xtick.major.size'] = 4
    plt.rcParams['ytick.major.size'] = 4
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['ytick.minor.size'] = 4
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['mathtext.fontset'] = 'dejavuserif' #'stix'
    plt.rcParams['font.family'] = 'serif' #'STIXGeneral'
    plt.rcParams['legend.facecolor'] = 'white'
    plt.rcParams['legend.fontsize'] = ticksize
    plt.rcParams['axes.formatter.limits'] = (-2,4)
    plt.rcParams['axes.linewidth'] = 1.75
    if grid:
        plt.rcParams['axes.grid'] = True
        plt.rcParams['grid.alpha'] = 0.6
    plt.rcParams['figure.figsize'] = (fxsize,fysize)
    plt.rcParams['figure.autolayout'] = True
    cm = plt.get_cmap('cividis')
    #cm = get_continuous_cmap(['#08232D', '#2D6C84', '#59982B', '#9FC126', '#DEEEA1'])
    #cm = mcolors.LinearSegmentedColormap.from_list('my_cmp', ['#08232D', '#2D6C84', '#59982B', '#9FC126', '#DEEEA1'], N=255)
    plt.rcParams['image.cmap'] = 'cividis'
    n_colours = 15
    lst = ['-', '--', '-.']#, ':']
    plt.rcParams['axes.prop_cycle'] = cycler(linestyle = n_colours//len(lst) * lst) + cycler(color = [cm(1.*i/n_colours) for i in range(n_colours)])
    
def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp