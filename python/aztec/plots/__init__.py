from aztec import *
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plotAstroSignal (axis, signal, ra, dec, range = None, **kwargs):
    if range is None:
        vmin,vmax = npy.percentile (signal, [0.1,99.9])
    else:
        vmin, vmax =range
    outImage = axis.imshow(npy.fliplr(signal.T), extent = [ra[-1],ra[0],dec[0],dec[-1]], origin = "lower", vmax = vmax, vmin=vmin, **kwargs)
    #outImage = axis.imshow(npy.fliplr(signal.T), extent = [ra.min(),ra.max(),dec[0],dec[-1]], origin = "lower", vmax = vmax, vmin=vmin, **kwargs)
    #outImage = axis.imshow(signal.T, extent = [ra[-1],ra[0],dec[0],dec[-1]], origin = "lower", vmax = vmax, vmin=vmin)
    
    return outImage

def plotMapInAxis (axis, aztecmap, range = None, units = units.mJy, cunits = units.arcmin, **kwargs):
    im= plotAstroSignal(axis, aztecmap.fSignal.to(units).value, aztecmap.RaCoords.to(cunits).value, aztecmap.DecCoords.to(cunits).value, range=range, **kwargs)
    axis.set_adjustable('box-forced')
    return im



def plotMapCountours(ax, amap, radii, levels):
  xx,yy = amap.generateCoordsMatrices()
  xx = xx.to("arcmin")
  yy = yy.to("arcmin")

  #
  wx = npy.where (npy.abs(amap.RaCoords.to("arcmin").value) <= radii)
  wy = npy.where (npy.abs(amap.DecCoords.to("arcmin").value) <= radii)
  XC = xx[wx[0][0]:wx[0][-1],wy[0][0]:wy[0][-1]]
  YC = yy[wx[0][0]:wx[0][-1],wy[0][0]:wy[0][-1]]
  pdatax = amap.fSignal.to("mJy").value[wx[0][0]:wx[0][-1],wy[0][0]:wy[0][-1]]

  cxr = ax.contourf(XC,YC,pdatax,levels=levels)
  #colorbar_inset(ax,cxr, nticks=6)
  ax.set_xlim(radii,-radii)
  ax.set_ylim(-radii, radii)
  
  return cxr


class CoordFormatter:
    def __init__(self,imap, ra, dec, units = ""):
        self.z = npy.fliplr(imap.T)
        self.x0 = ra[0]
        self.x1 = ra[-1]
        self.dy = ra[1]-ra[0]
        self.y0 = dec[0]
        self.y1 = dec[-1]
        self.dx = dec[1]-dec[0]
        self.units = units
        self.decs,self.ras = self.z.shape  #y,x
        print self.z.shape
        
    
    def format_coords(self,x,y):
        col = int(self.ras-(x-self.x0)/(self.x1-self.x0)*self.ras)
        row = int((y-self.y0)/(self.y1-self.y0)*self.decs)
        if col>=0 and col<self.ras and row>=0 and row<self.decs:
            z = self.z[row,col]
            return 'dra=%8.2f, ddec=%8.2f, value=%8.2f%s (%d,%d)' % (x, y, z, self.units,x,y)
        else:
            return 'x=%1.4f, y=%1.4f'%(col, row)



class MapPlotter:
    def __init__(self, fmap):
        self.fmap = fmap
        self.figure = None
        self.imageAxis = None
        self.range = None
        self.cbarAxis = None
        self.image = None
        self.lastSignal = None
        self.lastFiltered = None
        self.lastUnits = None
    
    def createFigureElements(self):
        self.figure = plt.figure(figsize=(5.5,5))
        self.imageAxis = plt.gca()
        lt1 = make_axes_locatable(self.imageAxis)
        self.cbarAxis = lt1.append_axes("right", size="5%", pad=0.3)
        
        
    def plot(self, variable = "signal", filtered = True, coordsUnits = units.arcsec):
        if self.figure is None:
            self.createFigureElements()
        else:
            try:
                self.figure.show()
            except:
                self.createFigureElements()
        
        signal, signal_units = self.selectSignal(variable,filtered)
        self.image = plotAstroSignal(self.imageAxis, signal, self.fmap.RaCoords.to(coordsUnits).value, \
                                     self.fmap.DecCoords.to(coordsUnits).value, range = self.range)
        formatter = CoordFormatter(signal, self.fmap.RaCoords.to(coordsUnits).value, self.fmap.DecCoords.to(coordsUnits).value,\
                                   units = signal_units)
        self.imageAxis.format_coord = formatter.format_coords
        self.imageAxis.set_xlabel(r"$\Delta \alpha$ [%s]"% coordsUnits.to_string())
        self.imageAxis.set_ylabel(r"$\Delta \delta$ [%s]"% coordsUnits.to_string())
        self.imageAxis.set_title("%s" % self.fmap.source)
        #Now create the colorbar
        self.figure.colorbar(self.image, cax=self.cbarAxis, orientation='vertical')
        self.cbarAxis.set_ylabel("%s" % signal_units)
        self.figure.show()
        self.lastFiltered=filtered
        self.lastSignal=variable
        self.lastUnits = coordsUnits
    
    def update(self):
        if not self.figure is None:
            self.plot(variable = self.lastSignal, filtered = self.lastFiltered, coordsUnits=self.lastUnits)
            self.figure.canvas.draw()
    def setRange (self, min,max):
        self.range =(min,max)
        self.update()

    def zoom (self, scale):
		try:
			scale = scale.to(self.lastUnits)
			scale = scale.value
		except AttributeError:
			pass
			
		
		self.imageAxis.set_xlim(scale,-scale)
		self.imageAxis.set_ylim(-scale, scale)
		self.update()

        
    def selectSignal (self,variable, filtered):
        oSignal = self.fmap.signal.to("mJy")
        sUnits = "Flux"
        
        if variable == "signal" and filtered:
            oSignal = self.fmap.fSignal.to("mJy")
        if variable == "weight":
            sUnits = "Weight"
            if filtered:
                oSignal = self.fmap.fWeight
            else:
                oSignal = self.fmap.weight
        if variable == "snr":
            sUnits = "SNR"
            if filtered:
                oSignal = self.fmap.fSnMap
            else:
                oSignal = self.fmap.snmap


        try:
			if sUnits == "Flux":
				sUnits +=  " [%s/beam]" %oSignal.unit.to_string()
			elif sUnits == "SNR":
				pass
			else:
				sUnits +=  " [%s]" %oSignal.unit.to_string()
			oSignal = oSignal.value
        except:
            pass

        return oSignal,sUnits
        
