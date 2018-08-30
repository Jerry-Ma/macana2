from lxml import etree
import xml.etree.ElementTree as ET
import os
import shutil
from aztec.macana.bstats import createBstatsFile


class apError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__ (self):
        return "Analysis Parameters Error: " + self.value

class apFileCreator:
    
    """This class allows to modify the content of an Analysis Parameters file (apfile). It is intended to automatize macana
       from python.
    """
    
    def __init__ (self, apTemplate):
        """Creates the apFileCreator
            Input:
                    apTemplate    An xml file with the Analysis Parameters you want to modiyf
        """
        self.template = apTemplate
        self.readTemplate()
        
    def readTemplate(self):
        parser = etree.XMLParser(remove_blank_text=True)
        self.apTree = etree.parse(self.template,parser=parser)
        self.analysisRoot = self.apTree.getroot()
        self.params = self.analysisRoot.find('parameters')
        self.coadd = self.analysisRoot.find('coaddition')
        self.noise = self.analysisRoot.find('noiseRealization')
        self.obs = self.analysisRoot.find('observations')
        self.simulate = self.analysisRoot.find('simulate')
        self.subtract = self.analysisRoot.find('subtract')
        self.filter = self.analysisRoot.find('wienerFilter')
        
    def changeParam(self, paramName, value):
        cElement = self.params.find(paramName)
        if cElement is None:
            raise apError('No such parameter on template file')
        cElement.text = str(value)
        
    def changeAttribute (self, paramName, attributeName, value):
        cElement = self.params.find(paramName)
        if cElement is None:
            raise apError('No such parameter on template file')
        try:
            cElement.attrib[attributeName]=value
        except Exception:
            raise apError('No such attribute on template file')
        
    
    def changeCoaddPath(self, newPath):
        path = self.coadd.find('mapPath')
        if not path is None:
            path.text = newPath + os.path.sep
        else:
            raise apError ('mapPath element not found in template file')
    
    def getCoaddMap (self):
        return self.coadd.find('mapFile').text
    
    def getCoaddMapPath(self):
        return self.coadd.find('mapPath').text
    
    def changeCoaddMap (self, newMapName):
        path = self.coadd.find('mapFile')
        if not path is None:
            path.text = newMapName
        else:
            raise apError ('mapPath element not found in template file')
        
    def changeSubParams (self, path, newMapName):
        if self.subtract is None: 
            self.subtract = etree.SubElement(self.analysisRoot, 'subtract')
            subPath = etree.SubElement(self.subtract, 'subPath')
            subFile = etree.SubElement (self.subtract, 'subFile')
        else:
            subPath = self.subtract.find('subPath')
            subFile = self.subtract.find('subFile')
        subPath.text = path
        subFile.text = newMapName
        
    def getSubParams(self):
        return self.subtract.find('subPath').text + os.path.sep + self.subtract.find('subFile').text
        
    def changeSimMap(self, newMapName):
        if not self.simulate is None:
            simMapE = self.simulate.find('simFile')
            simMapE.text = newMapName
        else:
            print "Parameter file has no simulate option, skipping"
    
    def changeSimParam(self, parameter, value):
        if not self.simulate is None:
            par = self.simulate.find (parameter)
            if not par is None:
                par.text = str(value)
            else:
                print "Not such simulation parameter, skipping"
        else:
            print "Parameter file has no simulate option, skipping"

    def changeNoisePath(self, newMapName):
        if not self.noise is None:
            noisePath = self.noise.find('noisePath')
            noisePath.text = newMapName
        else:
            print "Parameter file has no simulate option, skipping"
            
    def getNoisePath(self):
        if not self.noise is None:
            noisePath = self.noise.find('noisePath')
            return noisePath.text 
        else:
            print "Parameter file has no noise Path statement"

    def writeToFile(self, filename):
        f=open (filename, "w")
        f.write(etree.tostring(self.analysisRoot, pretty_print=True))
        f.close()
    
    def createSimulate (self, directory, simMap):
        if not self.simulate is None:
            self.analysisRoot.remove(self.simulate)
        self.simulate = ET.SubElement(self.analysisRoot, "simulate")
        fluxFactor = ET.SubElement(self.simulate, "fluxFactor")
        fluxFactor.text = "1"
        addSignal = ET.SubElement(self.simulate,"addSignal")  #  <addSignal>0</addSignal>
        addSignal.text= "0"
        atmFreq = ET.SubElement(self.simulate, "atmFreq")
        atmFreq.text = "0"
        atmFreq = ET.SubElement(self.simulate, "atmSeed")
        atmFreq.text = "-1"
        noiseChunk = ET.SubElement (self.simulate, "noiseChunk")
        noiseChunk.text = "0"
        simPath = ET.SubElement(self.simulate, "simPath")
        simPath.text = directory + os.path.sep
        simFile = ET.SubElement(self.simulate,"simFile")
        simFile.text = simMap
    
    def createSubtract(self, directory, subMap):
        if not self.subtract is None:
            self.analysisRoot.remove(self.subtract)
        self.subtract = ET.SubElement(self.analysisRoot, "subtract")
        subPath = ET.SubElement(self.subtract, "subPath")
        subPath.text = directory
        subFile = ET.SubElement(self.subtract,"subFile")  #  <addSignal>0</addSignal>
        subFile.text= subMap

    def setSingleObs (self, directory, ncfile, bstat_file = None, point_file = None):
        raw_data = directory +'/'
        noise_dir = directory + '/noise_maps/'
        map_dir = directory + '/reduced_maps/'
        self.changeCoaddPath(map_dir)
        self.changeNoisePath(noise_dir)
        self.makeObservations(raw_data, bstat_file=bstat_file, point_file = point_file, ncfile =ncfile)
        
    def setWorkingDir (self, directory, bstat_file = None, point_file = None):
        raw_data = directory + '/raw_data/'
        noise_dir = directory + '/noise_maps/'
        map_dir = directory + '/reduced_maps/'
        self.changeCoaddPath(map_dir)
        self.changeNoisePath(noise_dir)
        self.makeObservations(raw_data, bstat_file=bstat_file, point_file = point_file)

    def makeObservations (self, obsDir, bstat_file = None, point_file = None, ncfile = None):
        
        if ncfile is None:
            import fnmatch
            fileList=[]
            files = os.listdir(obsDir)
            for ifile in files:
                if fnmatch.fnmatch(ifile, "*.nc"):
                    fileList.append(os.path.join(obsDir,ifile))
            fileList.sort()
        else:
            fileList = [ncfile]
        nfiles = len(fileList)
        if nfiles > 0:
            if self.obs is None:
                self.obs = etree.Element("observations")
            else:
                self.analysisRoot.remove(self.obs)
                self.obs = etree.Element("observations")
            
            if point_file is None:
                pointingCor = None
            else:
                from aztec.macana.pointing import LMTPointingTextFile
                pointingCor = LMTPointingTextFile(point_file)
            
            fileList.sort()
            rawDir = etree.Element("rawDataPath")
            rawDir.text = obsDir
            bsPath = etree.Element("bsPath")
            bsPath.text = obsDir
            mapPath = etree.Element("mapPath")
            mapPath.text = self.coadd.find('mapPath').text
            nFiles = etree.Element("nFiles")
            nFiles.text = str(nfiles) 
            
            self.obs.append(rawDir)
            self.obs.append(bsPath)
            self.obs.append(mapPath)
            self.obs.append(nFiles)
            curfile = 0
            
            print "Adding %d files to parameter file" % nfiles
            for ifile in fileList:
                cFile = etree.Element ("f"+str(curfile))
                strFile = etree.Element("fileName")
                filename = os.path.basename(ifile)
                strFile.text = filename
                mapFile = etree.Element("mapName")
                ofile = filename.replace(".nc", "_map.nc")
                mapFile.text=ofile
                bstat= etree.Element("bsName")
                if bstat_file is None:
                    ind_bstat = os.path.basename(ifile.replace('.nc', '.bstats'))
                    if not os.path.exists(os.path.join(obsDir, ind_bstat)):
                        bfile = createBstatsFile(ifile)
                        if not os.path.exists(obsDir + bfile):
                            shutil.copy2(os.path.join (os.path.expanduser(os.environ['AZTEC_MACANA_PATH']),"default_bstats",bfile), obsDir+bfile)
                    else:
                        bfile = ind_bstat
                else:
                    bfile = bstat_file
                bstat.text = bfile
                if not pointingCor is None:
                    try:
                        from aztec.util import getObsNum
                        obsnum = getObsNum(ifile)
                        azOff = -1.0*pointingCor.az[obsnum]
                        elOff = -1.0*pointingCor.el[obsnum]
                    except KeyError:
                        print "Warning: No pointing correction for file %s" % ifile
                        azOff = 0.0;
                        elOff = 0.0;
                else:
                    azOff = 0.0;
                    elOff = 0.0;
                boff0 = etree.Element("bsOffset_0")
                boff0.text = "%10.4f" % azOff
                boff1 = etree.Element("bsOffset_1")
                boff1.text =  "%10.4f" % elOff  
                cFile.append(strFile)
                cFile.append(bstat)
                cFile.append(mapFile)
                cFile.append(boff0)
                cFile.append(boff1)
                self.obs.append(cFile)
                curfile+=1

            self.analysisRoot.append(self.obs)
    
