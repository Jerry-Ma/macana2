from aztec import *
from scipy.io.netcdf import netcdf_file as NetCDFFile
from scipy.io.idl import readsav


class Map:
    """Base Class for any AzTEC Map Structure"""
    
    
    
    def __init__ (self,signal, ra, dec, raCenter = 0, decCenter = 0, source="Source"):
        self.signal = signal
        self.source = source
        self.RaCoords = ra
        self.DecCoords = dec
        self.RaCenter = raCenter 
        self.DecCenter = decCenter
    
    def generateRGrid(self, a=1.0,b=1.0,angle=0.0, x0=0., y0=0.):
        """Generate a matrix object with the distance relative to the map center"""
        xx,yy = self.generateCoordsMatrices(x0=x0,y0=y0)
        xx/=a
        yy/=b
        if angle != 0:
            print angle
            rangle = npy.radians(angle)
            x = xx*npy.cos(rangle) - yy*npy.sin(rangle)
            y = xx*npy.sin(rangle) + yy*npy.cos(rangle)
            xx = x
            yy = y
        return npy.sqrt(xx**2 + yy**2)
    def generateAngleGrid(self, a=1.0,b=1.0,angle=0.0, x0=0., y0=0.):
        """Generate a matrix object with the distance relative to the map center"""
        xx,yy = self.generateCoordsMatrices(x0=x0,y0=y0)
        if angle != 0:
            rangle = npy.radians(angle)
            x = xx*npy.cos(rangle) - yy*npy.sin(rangle)
            y = xx*npy.sin(rangle) + yy*npy.cos(rangle)
            xx = x
            yy = y
        return npy.arctan2(yy.value,xx.value)
    
    def generateCoordsMatrices(self, x0=0.0, y0=0.0):
        xx,yy = npy.meshgrid((self.RaCoords-x0),(self.DecCoords-y0))
        return xx.T, yy.T

class AztecMap(Map):
    
    def __init__ (self,ncFileName, source=None):
        if os.path.exists(ncFileName):
            self.filename = ncFileName
            self.filtered = False
            if not source is None:
                self.source = source
            try:
                ncFile = NetCDFFile(ncFileName)
                self.loadFromNcFile(ncFile)
                ncFile.close()
            except RuntimeError:
                #try:
                savFile = readsav (ncFileName)
                self.loadFromIDLSav(savFile)
                #except Exception:
                    #raise  Exception ("AzTEC Map Error. Unknown format")
#             
            self.fixWeightMap()
            self.calculateSNMap()
            
        else:
            raise Exception("AztecMap Error. No such file or directory")
    
    def loadFromNcFile(self, ncFile):
            
            """Create a AzTEC Map object from a nc file. This is aztec c++ pipeline output"""    
            try:
                self.source = getattr(ncFile,"source").replace("\"","")
            except AttributeError:
                self.source = ""
                pass    
            try:
                self.timeChunk = ncFile.timeChunk
            except AttributeError:
                self.timeChunk = -1
            try:
                self.signal = npy.array(ncFile.variables['signal'][:])*units.Jy
            except KeyError:
                self.signal = npy.array(ncFile.variables['noise'][:])*units.Jy
            try:
                self.weight = npy.array(ncFile.variables['weight'][:])* (1.0/units.Jy)**2.0
            except KeyError:
                self.weight = None

            try:
                self.kernel = npy.array(ncFile.variables['kernel'][:])
            except KeyError:
                self.kernel = None
            try:
                self.tSignal = npy.array(ncFile.variables['tSignal'][:])*units.Jy
            except KeyError:
                self.tSignal = None
            try:
                self.tWeight = npy.array(ncFile.variables['tWeight'][:])* (1.0/units.Jy)**2.0
            except KeyError:
                self.tWeight = None
            try:
                self.tKernel = npy.array(ncFile.variables['tKernel'][:])
            except KeyError:
                self.tKernel = None
            self.RaCoords = npy.array(ncFile.variables['rowCoordsPhys'][:])*180.0/(npy.pi)*units.deg
            self.DecCoords = npy.array(ncFile.variables['colCoordsPhys'][:])*180.0/(npy.pi)*units.deg
            self.AbsRaCoords = npy.array(ncFile.variables['xCoordsAbs'][:])*180.0/(npy.pi)*units.deg
            self.AbsDecCoords = npy.array(ncFile.variables['yCoordsAbs'][:])*180.0/(npy.pi)*units.deg
            try:
                self.kernel = npy.array(ncFile.variables['kernel'][:])
            except KeyError:
                self.kernel = None
            try:
                self.SourceRa = getattr(ncFile, "MasterGrid[0]")*180.0/npy.pi*units.deg
                self.SourceDec = getattr(ncFile, "MasterGrid[1]")*180.0/npy.pi*units.deg
            except AttributeError:
                self.SourceRa = getattr(ncFile, "RaCenter")*180.0/npy.pi*units.deg
                self.SourceDec = getattr(ncFile, "DecCenter")*180.0/npy.pi*units.deg
            
            self.RaCenter = 0.0* units.deg
            self.DecCenter = 0.0* units.deg
            try:
                self.fSignal = npy.array(ncFile.variables['filteredSignal'][:])*units.Jy
                self.fWeight = npy.array(ncFile.variables['filteredWeight'][:])* (1.0/units.Jy)**2.0
                self.fKernel = npy.array(ncFile.variables['filteredKernel'][:])
                self.filtered = True
            except KeyError:
                self.fSignal = None
                self.fWeight = None
                self.filtered = False
            try: 
                self.psd2d = npy.array(ncFile.variables['psd_2d'][:])
                self.psdFreq2d = npy.array(ncFile.variables['psdFreq_2d'][:])
            except KeyError:
                self.psd2d = None
                self.psdFreq2d = None
            
            try:
                self.atmMap = npy.array(ncFile.variables['atmMap'][:])
            except KeyError:
                self.atmMap = None
            
            self.sourceArray = []
            try:
                self.nSources = len(ncFile.dimensions['nSources'])
            except KeyError:
                self.nSources = 0
	    except TypeError:
		self.nSources = ncFile.dimensions['nSources']
            
            try:
                cflux = ncFile.variables['signal'].getncattr ("amplitude")
                cflux_err =  ncFile.variables['signal'].getncattr ("amplitude_err")
                cFWHM_x = ncFile.variables['signal'].getncattr ("FWHM_x")
                cFWHM_x_err = ncFile.variables['signal'].getncattr ("FWHM_x_err")
                cFWHM_y = ncFile.variables['signal'].getncattr ("FWHM_y")
                cFWHM_y_err = ncFile.variables['signal'].getncattr ("FWHM_y_err")
                cOffset_x = ncFile.variables['signal'].getncattr ("offset_x")
                cOffset_x_err = ncFile.variables['signal'].getncattr ("offset_x_err")
                cOffset_y = ncFile.variables['signal'].getncattr ("offset_y")
                cOffset_y_err = ncFile.variables['signal'].getncattr ("offset_y_err")
                self.coaddFit={'flux':cflux, 'flux_err':cflux_err, 'FWHM_x':cFWHM_x,'FWHM_y':cFWHM_y, \
                               'FWHM_x_err':cFWHM_x_err, 'FWHM_y_err':cFWHM_y_err,\
                               'Offset_x': cOffset_x, 'Offset_y': cOffset_y,\
                               'Offset_x_err': cOffset_x_err, 'Offset_y_err': cOffset_y_err}
            except KeyError:
                self.coaddFit=None
            except AttributeError:
                self.coaddFit=None
                
            if self.nSources > 0:
                for i in range (self.nSources):
                    svar = "pointSource_%d" % i
                    centerFlux = ncFile.variables[svar].centerFlux
                    dcFlux = ncFile.variables[svar].dc_offset
                    flux = ncFile.variables[svar].amplitude
                    #xpos = ncFile.variables[svar].centerXPos
                    #ypos = ncFile.variables[svar].centerYPos
                    xpos = ncFile.variables[svar].raPhysCentroid
                    ypos = ncFile.variables[svar].decPhysCentroid
                    rapos = ncFile.variables[svar].offset_x
                    decpos = ncFile.variables[svar].offset_y
                    sdict = {'flux':flux, 'cflux':centerFlux, 'dcOffset':dcFlux, 'xpos':xpos, 'ypos':ypos, 'rapos': rapos, 'decpos': decpos}
                    self.sourceArray.append(sdict)

    
    def loadFromIDLSav(self,savHandler):
            self.filtered = True
            self.signal = npy.array(savHandler.coadded_signal_map, copy = True).T*units.Jy
            self.fSignal = self.signal
            self.weight = npy.array(savHandler.coadded_weight_map, copy = True).T* (1.0/units.Jy)**2.0
            self.fWeight = self.weight
            self.SourceRa = savHandler.az_coord_center
            self.SourceDec = savHandler.el_coord_center
            self.RaCoords = npy.array(savHandler.phys_ra_map, copy = True)
            self.RaCoords = (self.RaCoords-self.SourceRa)*units.deg
            self.DecCoords = npy.array(savHandler.phys_dec_map, copy = True)
            self.DecCoords = (self.DecCoords-self.SourceDec)*units.deg
            self.AbsRaCoords = npy.array(savHandler.abs_ra_map,copy = True).T
            self.AbsDecCoords = npy.array(savHandler.abs_dec_map,copy = True).T
            self.kernel = npy.array(savHandler.kernel_mean).T
            self.fKernel = self.kernel
            self.valid = True
            self.source = ''
            self.coaddFit = None
            #self.fixKernel()

    def isFiltered(self):
        return self.filtered
    
    def fixWeightMap(self):
        self.zeroNegativeValues(self.weight)
        if self.filtered:
            self.zeroNegativeValues(self.fWeight)
            
    def zeroNegativeValues (self, var):
        
        if var is None:
            return
        if not isinstance( var, ( npy.ndarray) ):
            return
        var = npy.where(var < 0.0, 0.0, var)
        
    def findBounds(self,coords = False):
        #Return the pixel bounds for the map
        ix = 0
        ex = -1
        if self.filtered:
            signal = self.fWeight.value
        else:
            signal = self.weight.value
        ix = signal.shape[1]/2
        ex = signal.shape[1]/2
        while (sum(abs(signal[:,ix])) != 0):
            ix-=1
        while (sum(abs(signal[:,ex])) != 0):
            ex+=1
        iy = signal.shape[0]/2
        ey = signal.shape[0]/2    
        while (sum(abs(signal[iy,:])) != 0):
            iy-=1
        while (sum(abs(signal[ey,:])) != 0):
            ey+=1

        if coords:
            return [self.RaCoords[ix+1], self.RaCoords[ex-1], self.DecCoords[iy+1], self.DecCoords[ey-1]]
        else:
            return [ix+1,ex-1,iy+1,ey-1]
    
    
    def calculateSNMap(self):
        self.snmap =  self.signal * npy.sqrt(self.weight)
        if self.filtered:
            self.fSnMap = self.fSignal * npy.sqrt(self.fWeight)
    
    def getCovCutIndex (self, wcut, filtered = True):
        if self.filtered and filtered:
            wmap = self.fWeight
        else:
            wmap = self.weight
        
        covcut = npy.percentile(wmap[npy.where(wmap>0.0)], wcut)
        wremove = npy.where(wmap.value < covcut)
        winclude = npy.where(wmap.value >=covcut)
        
        return wremove,winclude
    
    def wcut (self,wcut):
        """This removes all pixel from signal,weight and SN maps that are below the wcut percentile value on the weight map
            Input:
                    wcut    Threshold for cutting, all the pixels below this percentile value from 
                            the weight map will be set to 0.0
        """
        if self.filtered:
            wmap = self.fWeight
        else:
            wmap = self.weight
        
        w,_=self.getCovCutIndex(wcut, filtered=self.filtered)
        
        if len (w) > 0:
            self.signal[w]=0.0
            self.weight[w]=0.0
            if not isinstance(self.kernel, int):
                self.kernel[w]=0.0
            if self.filtered:
                self.fSignal[w]=0.0
                self.fWeight[w]=0.0
                self.fKernel[w]=0.0
            self.calculateSNMap()
            
    def snCut (self,sncut):
        """This removes all pixel from signal and SN maps that are below the sncut value on the SN map
            Input:
                    sncut    Threshold for cutting, all the pixels below this value from 
                            the SN map will be se to 0.0
        """
        if self.filtered:
            snmap = self.fSnMap
        else:
            snmap = self.snmap
        w = npy.where(snmap < sncut)
        if len(w) > 0:
            self.signal[w]=0.0
            if self.filtered:
                self.fSignal = 0.0
            self.calculateSNMap()
    
    def writeToNc(self,ncFileName):
        ncFile = NetCDFFile(ncFileName, 'w')
        ncFile.createDimension("nrows", self.RaCoords.shape[0])
        ncFile.createDimension("ncols", self.DecCoords.shape[0])
        dims = ("nrows","ncols")
        drow = ("nrows",)
        dcol = ("ncols",)
        ncFile.createVariable("signal", "d", dims)
        ncFile.createVariable("weight", "d", dims)
        ncFile.createVariable("kernel", "d", dims)
        ncFile.createVariable("rowCoordsPhys", "d", drow)
        ncFile.createVariable("colCoordsPhys", "d", dcol)
        ncFile.createVariable("xCoordsAbs", "d", dims)
        ncFile.createVariable("yCoordsAbs", "d", dims)
        ncFile.createVariable("filteredSignal", "d", dims)
        ncFile.createVariable("filteredWeight", "d", dims)
        ncFile.createVariable("filteredKernel", "d", dims)
        setattr(ncFile, "source", "%s"%self.source)
        setattr(ncFile, "MasterGrid[0]", self.SourceRa/180.0*(npy.pi))
        setattr(ncFile, "MasterGrid[1]", self.SourceDec/180.0*(npy.pi))
        ncFile.variables['signal'][:] = self.signal.value
        ncFile.variables['weight'][:] = self.weight.value
        ncFile.variables['rowCoordsPhys'][:] = self.RaCoords /180.0*(npy.pi)
        ncFile.variables['colCoordsPhys'][:] = self.DecCoords /180.0*(npy.pi) 
        ncFile.variables['kernel'][:] = self.kernel
        ncFile.variables['xCoordsAbs'][:] = self.AbsRaCoords.to("rad").value
        ncFile.variables['yCoordsAbs'][:] = self.AbsDecCoords.to("rad").value
        if self.filtered:
            ncFile.variables['filteredSignal'][:] = self.fSignal.value
            ncFile.variables['filteredWeight'][:] = self.fWeight.value
            ncFile.variables['filteredKernel'][:] = self.fKernel
        else:
            ncFile.variables['filteredSignal'][:] = self.signal.value
            ncFile.variables['filteredWeight'][:] = self.weight.value
            ncFile.variables['filteredKernel'][:] = self.kernel
        ncFile.sync()    
        ncFile.close()
