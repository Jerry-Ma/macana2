from aztec import *
from numpy import fft as FFT


class KernelError(Exception):
    def __init__(self):
        super(KernelError, self).__init__("Kernel and Signal shape miss match")

class Filter:
    
    def __init__(self, kernel):
        self.kernel = kernel
        self.fkernel = FFT.fft2(kernel)
        
    def filter (self,signal):
        if signal.shape[0] != self.kernel.shape[0] or signal.shape[1] != self.kernel.shape[1]:
            raise KernelError()
        signalTmp = signal.copy()
        fSignal = FFT.fft2(signal)
        fSignal*=self.fkernel
        return FFT.ifftshift(FFT.ifft2(fSignal)).real

class GaussianFilter(Filter):
    
    def __init__(self,x,y,pars,normalize=1.0):
        
        from aztec.gaussfit import gaussian2d
        
        self.xx, self.yy = npy.meshgrid(x,y)
        self.xx = self.xx.transpose()
        self.yy = self.yy.transpose()
        kernel = gaussian2d((self.xx,self.yy), *pars)/normalize
        Filter.__init__(self, kernel)
        
class MapGaussianFilter(GaussianFilter):
    def __init__(self,ufmap,beamsize, normalize = "pix", replace = True):
        self.ufmap = ufmap
        self.beamsize = beamsize
        if isinstance(self.beamsize, (int,long,float)):
            self.beamsize *=units.deg
        self.sigma = self.beamsize*beam2sigma
        self.pixSize = npy.abs(ufmap.RaCoords[0]-ufmap.RaCoords[1])
        self.setNormalizationValue(normalize)
        self.replace = replace
        pars = [0.0, 1.0, self.getCenterValue(ufmap.RaCoords).value, \
                self.getCenterValue(ufmap.DecCoords).value, self.sigma.to(ufmap.RaCoords.unit).value, \
                self.sigma.to(ufmap.DecCoords.unit).value, 0.0]
        
        print pars
        GaussianFilter.__init__(self, ufmap.RaCoords, ufmap.DecCoords, pars, normalize = self.normvalue)
    
    def doFilter(self):
        from scipy.ndimage.filters import gaussian_filter
        import warnings
        signalMap = self.ufmap.signal
        kernelMap = self.ufmap.kernel
        wMapMask = npy.where (self.ufmap.weight>0.0, True, False)
        warnings.simplefilter("ignore")
        noiseMap = npy.where (wMapMask, 1.0/npy.sqrt(self.ufmap.weight), 0.0)
        #noiseMap /= (self.beamsize.to (self.pixSize.unit)/self.pixSize)      
        fSignal = self.filter(signalMap)
        fNoise = self.filter(noiseMap)
        fKernel = self.filter(kernelMap)
        #Reduce the filtered noise by the beamsize ratio, this assumes gaussian noise only
        fWeight = npy.where (wMapMask, 1.0/fNoise**2.0, 0.0)

        if not self.replace:
            return fSignal, fWeight, fKernel
        else:
            self.ufmap.fSignal = fSignal * self.ufmap.signal.unit
            self.ufmap.fWeight = fWeight * self.ufmap.weight.unit
            self.ufmap.fKernel = fKernel
            self.ufmap.filtered = True
            self.ufmap.calculateSNMap()
        
    def setNormalizationValue (self,normalize):
        self.normvalue = 1.0
        if normalize == "pix":
            sigma = (self.beamsize/self.pixSize)*beam2sigma
            self.normvalue = 2.0*npy.pi*(sigma**2.0)
            self.normvalue = self.normvalue.decompose().value
        elif normalize == "area":
            sigma = self.beamsize
            self.normvalue = 2.0*npy.pi*(sigma**2.0)
            self.normvalue = self.normvalue.decompose().value
    
    def getCenterValue(self,coordinate):
        varSize = coordinate.shape[0]
        center = int(varSize/2)
        isOdd = varSize % 2
        
        if isOdd == 1:
            return coordinate[center]
        else:
            return 0.5*(coordinate[center]+coordinate[center+1]) 