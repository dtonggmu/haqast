#!/data/aqf/barryb/anaconda2/bin/python
from numpy import *
import matplotlib.pyplot as plt
import seaborn
import pytz
import xarray as xr
from datetime import datetime
import sys
from numpy import arange,array

def readfile(file):
    from netCDF4 import MFDataset
    return MFDataset(file)

def make_spatial_plot(cmaqvar, x, y, date, m, dpi=None, savename='', levs=arange(10,110,10), cmap='YlGnBu'):
    from numpy import arange
    fig = plt.figure(figsize=(18, 10), frameon=False)
    # define map and draw boundries
    m.drawstates()
    m.drawcoastlines(linewidth=.3)
    m.drawcountries()
    plt.axis('off')
    ncolors = len(levs)
    c, cmap = colorbar_index(ncolors, cmap,levs)
    m.pcolormesh(x, y, cmaqvar, vmin=min(levs), vmax=max(levs), cmap=cmap)
    titstring = date
    plt.title(titstring)

    plt.tight_layout()
    if savename != '':
        plt.savefig(savename + date.strftime('%Y%m%d_%H.jpg'), dpi=dpi)
        plt.close()
    return c
 
def colorbar_index(ncolors, cmap, levels):
    import matplotlib.cm as cm
    import numpy as np
    cmap = cmap_discretize(cmap, ncolors)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors + 0.5)
    colorbar = plt.colorbar(mappable,format='%1.2g')
    colorbar.set_ticks(np.linspace(0,ncolors,ncolors))
    colorbar.set_ticklabels(levels)
    
    return colorbar, cmap

  
def cmap_discretize(cmap, N):
    """
    Return a discrete colormap from the continuous colormap cmap.
    cmap: colormap instance, eg. cm.jet. 
    N: number of colors.
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
    import numpy as np

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in xrange(N + 1)]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)

def load_conus_basemap(grdobj):
    from mpl_toolkits.basemap import Basemap
    latitude = grid.variables['LAT'][:][0, 0, :, :].squeeze()
    longitude = grid.variables['LON'][:][0, 0, :, :].squeeze()
    lat1 = grid.P_ALP
    lat2 = grid.P_BET
    lon1 = grid.P_GAM
    lon0 = grid.XCENT
    lat0 = grid.YCENT
    m = Basemap(projection='laea', resolution='h', lat_1=lat1, lat_2=lat2, lat_0=lat0, lon_0=lon0, lon_1=lon1,
                    llcrnrlat=latitude[0, 0], urcrnrlat=latitude[-1, -1], llcrnrlon=longitude[0, 0],
                    urcrnrlon=longitude[-1, -1], rsphere=6371200.,
                    area_thresh=50.)
    x,y = m(longitude,latitude)
    return m,x,y
  
def get_dates(concobj):
    from datetime import datetime
    from pandas import DataFrame
       
    from numpy import concatenate, arange,array
    tflag1 = array(concobj.variables['TFLAG'][:, 0, 0], dtype='|S7')
    tflag2 = array(concobj.variables['TFLAG'][:, 1, 1] / 10000, dtype='|S6')
    date = []
    for i, j in zip(tflag1, tflag2):
        date.append(datetime.strptime(i + j, '%Y%j%H'))
    dates = array(date)
    r = DataFrame(dates, columns=['dates'])
    rr = r.drop_duplicates(keep='last')
    indexdates = rr.index.values
    utc = pytz.timezone('UTC')
    est = pytz.timezone('US/Eastern')
    dates = [ utc.localize(i) for i in dates]
    est_dates = array([i.astimezone(est) for i in dates])
    return est_dates[indexdates],indexdates
  
def main(concfile,gridfile,param):
    from netCDF4 import Dataset,MFDataset
    from mpl_toolkits.basemap import Basemap
    from numpy import unique,array
    
    concfiles = ['temp/20161114/aqm.t12z.aconc.ncf']
    grid = Dataset('MAY2014/aqm.t12z.grdcro2d.ncf')
    
    concobj = MFDataset(concfiles)
    print 'opened file'
    d,index = get_dates(concobj)
    print 'loading basemap'
    latitude = grid.variables['LAT'][:][0, 0, :, :].squeeze()
    longitude = grid.variables['LON'][:][0, 0, :, :].squeeze()
    lat1 = grid.P_ALP
    lat2 = grid.P_BET
    lon1 = grid.P_GAM
    lon0 = grid.XCENT
    lat0 = grid.YCENT
    m = Basemap(projection='laea', resolution='h', lat_1=lat1, lat_2=lat2, lat_0=lat0, lon_0=lon0, lon_1=lon1,
                    llcrnrlat=latitude[0, 0], urcrnrlat=latitude[-1, -1], llcrnrlon=longitude[0, 0],
                    urcrnrlon=longitude[-1, -1], rsphere=6371200.,
                    area_thresh=50.)
    x,y = m(longitude,latitude)

    
    #get ozone
    o3 = concobj.variables['O3'][index,0,:,:].squeeze() * 1000.
    #to make an image for each time loop through
    for i,j in enumerate(d):
      date = j.strftime('%m/%d/%Y %H')
      make_spatial_plot(o3[i,:,:], x, y, date, m, levs=arange(10,110,10), cmap='RdBu_r')
      plt.savefig(j.strftime('%Y%m%d%H_o3.jpg'),dpi=100)
      print ' plot made for: ' + date
      plt.close()
    
    
if __name__ == '__main__':
  main(sys.argv[1],sys.argv[2],sys.argv[3])
