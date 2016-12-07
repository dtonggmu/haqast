import matplotlib.colors as mcolors
from matplotlib import cm
from numpy import vstack,arange

def o3cmap():
    # This function returns the colormap and bins for the ozone spatial plots
    #this is designed to have a vmin =0 and vmax = 140
    #return cmap,bins
    colors1 = cm.viridis(linspace(0,1,128))
    colors2 = cm.OrRd(linspace(.2,1,128))
    colors = vstack((colors1,colors2))
    return mcolors.LinearSegmentedColormap.from_list('o3cmap',colors),arange(0,140.5,.5)

def pm25cmap():
    # This function returns the colormap and bins for the PM spatial plots
    #this is designed to have a vmin =0 and vmax = 140
    #return cmap,bins
    colors1 = cm.viridis(linspace(0,1,128))
    colors2 = cm.OrRd(linspace(.2,1,128))
    colors = vstack((colors1,colors2))
    return mcolors.LinearSegmentedColormap.from_list('pm25cmap',colors),arange(0,70.2,.2)

def noxcmap():
    # This function returns the colormap and bins for the NO2/NO/NOx spatial plots
    #this is designed to have a vmin =0 and vmax = 140
    #return cmap,bins
    colors1 = cm.viridis(linspace(0,1,128))
    colors2 = cm.plasma_r(linspace(.042,.75,128))
    colors = vstack((colors1,colors2))
    return mcolors.LinearSegmentedColormap.from_list('noxcmap',colors),arange(0,70.2,.2)

def so2cmap():
    # This function returns the colormap and bins for the NO2/NO/NOx spatial plots
    #this is designed to have a vmin =0 and vmax = 140
    #return cmap,bins
    colors1 = cm.viridis(linspace(0,1,128))
    colors2 = cm.plasma_r(linspace(.042,.75,128))
    colors = vstack((colors1,colors2))
    return mcolors.LinearSegmentedColormap.from_list('noxcmap',colors),arange(0,75.2,.2)

def pm10cmap():
    # This function returns the colormap and bins for the NO2/NO/NOx spatial plots
    #this is designed to have a vmin =0 and vmax = 140
    #return cmap,bins
    colors1 = cm.viridis(linspace(0,1,128))
    colors2 = cm.plasma_r(linspace(.042,.75,128))
    colors = vstack((colors1,colors2))
    return mcolors.LinearSegmentedColormap.from_list('noxcmap',colors),arange(0,150.5,.5)
