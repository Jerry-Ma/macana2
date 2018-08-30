from aztec import *

def getAstePointing (ncfile):
    from aztec.util import getAztecSeason, getJulDatefromNc
    from scipy.io.idl import readsav
    from scipy.interpolate import interp1d
    
    
    
    jdate = getJulDatefromNc(ncfile)
    season =getAztecSeason(jdate) 
    
    idlpath = os.getenv("AZTEC_IDL_PATH")
    boresight = os.path.join(idlpath, "parameters_%s" %season, "boresight_info.sav")
    
    dataDic = readsav(boresight)
    jdates = dataDic['boresight_info'][0][0]
    lst = dataDic['boresight_info'][0][1]
    azoff = dataDic['boresight_info'][0][2]
    eloff = dataDic['boresight_info'][0][3]
    
    wdate = npy.where(npy.logical_and(jdates >=jdate -0.5,jdates < jdate+0.5))
    
    jdates = jdates[wdate]
    lst = lst[wdate]
    azoff = azoff[wdate]
    eloff = eloff[wdate]
    
    if jdates.size == 0:
        raise Exception("No data to interpol")
    
    mylst = astro.time.Time(jdate, scale = "utc", format ="jd").sidereal_time("apparent", ASTE_position.longitude).value
     
    wdate = npy.where (npy.logical_and(lst >= mylst-0.5, lst < mylst+0.5))
     
    jdates = jdates[wdate]
    lst = lst[wdate]
    azoff = azoff[wdate]
    eloff = eloff[wdate]
     
    if jdates.size == 0:
        raise Exception("No data to interpol")
    
    if jdates.size == 1:
        return azoff[0],eloff[0]
    
    az = npy.interp (mylst, lst, azoff)
    el = npy.interp (mylst, lst, eloff)
    
    return az,el

def reduceSingleLMTPointing (ncfile, apFile = 'apLMTP.xml', addRaDec = False):
    from aztec.macana.aptools import apFileCreator
    from aztec.macana import runMacana
    from aztec.util import getObsNum, getJulDatefromNc, getOffsetFromMap
    
    fileDir = os.path.dirname(os.path.abspath(ncfile))
    ap = apFileCreator(apFile)
    ap.setSingleObs (fileDir, ncfile)
    apFileName = "ap%s.xml" % os.path.basename(ncfile).replace(".nc","")
    ap.writeToFile(apFileName)
    obsnum = getObsNum(ncfile)
    jdate = getJulDatefromNc(ncfile)
    dateid = npy.floor(jdate-0.5)
    runMacana (apFileName)
    outfile = os.path.join(ap.getCoaddMapPath(), ap.getCoaddMap())
    azoffset,eloffset = getOffsetFromMap(outfile)
    if addRaDec:
       ap.changeParam("azelMap",0)
       ap.observations.find("f0").find("bsOffset_0").text = str(-1.0*azOffset)
       ap.observations.find("f0").find("bsOffset_1").text = str(-1.0*elOffset)
       ap.writeToFile(apFileName)
       runMacana (apFileName)
       raoffset, decoffset = getOffsetFromMap(outfile)
       return dateid, obsnum, jdate, azoffset, eloffset, raoffset, decoffset
    
    return dateid, obsnum, jdate, azoffset, eloffset


def createPointingFile(pointingDir, apFile = 'apLMTP.xml', outfile = 'pointing_info.pys', addRaDec=False):
    import fnmatch
    from pickle import Pickler
    
    files = os.listdir(pointingDir)
    ncFiles = []
    for ifile in files:
        if fnmatch.fnmatch(ifile,"*.nc"):
            ncFiles.append(os.path.join(pointingDir,ifile))
    ncFiles.sort()
    pointingDict = {}
    for inc in ncFiles:
        filep = reduceSingleLMTPointing(os.path.abspath(inc),apFile=apFile, addRaDec=addRaDec)
        if not filep[0] in pointingDict.keys():
            pointingDict[filep[0]] = npy.array(filep[1:])
        else:
            pointingDict[filep[0]] = npy.vstack((pointingDict[filep[0]], filep[1:]))
    
    fhand = open(outfile,'w')
    phand = Pickler(fhand)
    phand.dump(pointingDict)
    fhand.close()


def getLMTPointing(rawDir, pointingDir, outfile = 'interpolated_pointing.txt', overwrite=False, apFile = 'apLMTP.xml', reduce=True):
    from pickle import Unpickler
    from scipy.interpolate import interp1d
    import fnmatch
    from aztec.util import getObsNum, getJulDatefromNc
    
    pinfoName = os.path.join(os.path.abspath(pointingDir),"pointing_info.pys")
    if reduce:
        createPointingFile(pointingDir, outfile=pinfoName, apFile=apFile)
    fhand = open(pinfoName,'r')
    phand = Unpickler(fhand)
    pDict = phand.load()
    fhand.close()
    #Create Interpolators
    intDict = {}
    for idate in pDict.keys():
        if len(pDict[idate].shape) == 1:
            pDict[idate] = pDict[idate].reshape([1,4])
	try:
		try:
		    print pDict[idate].shape
		    azInterp = interp1d(pDict[idate][:,1],pDict[idate][:,2], type ='quadratic')
		except TypeError:
		     azInterp = interp1d(pDict[idate][:,1],pDict[idate][:,2])
		try:
		    elInterp = interp1d(pDict[idate][:,1],pDict[idate][:,3], type ='quadratic')
		except TypeError:
		    elInterp = interp1d(pDict[idate][:,1],pDict[idate][:,3])
		jdmax = pDict[idate][:,1].max()
		jdmin = pDict[idate][:,1].min()
		intDict[idate]={'azInterp':azInterp, 'elInterp':elInterp, 'jdmax':jdmax, 'jdmin':jdmin, 'valid':True}
	except ValueError:
		jdmax = pDict[idate][:,1].max()
		jdmin = pDict[idate][:,1].min()
		intDict[idate]={'azInterp':pDict[idate][0,2], 'elInterp':pDict[idate][0,3], 'jdmax':jdmax, 'jdmin':jdmin, 'valid':False}
    
    files = os.listdir(rawDir)
    ncFiles = []
    for ifile in files:
        if fnmatch.fnmatch(ifile,"*.nc"):
            ncFiles.append(os.path.join(rawDir,ifile))
    ncFiles.sort()
    
    ofile = open (outfile,'w')
    for ifile in ncFiles:
        obsnum = getObsNum(ifile)
        jdate = getJulDatefromNc(ifile)
        jdateid = npy.floor(jdate-0.5)
        if not jdateid in intDict.keys():
            raise Exception ("No pointing information for obsnum %d" % obsnum)
        if jdate > intDict[jdateid]['jdmax']:
            jdate = intDict[jdateid]['jdmax']
        if jdate < intDict[jdateid]['jdmin']:
            jdate = intDict[jdateid]['jdmin']
	if intDict[jdateid]['valid']:
            azOffset = intDict[jdateid]['azInterp'](jdate)
            elOffset = intDict[jdateid]['elInterp'](jdate)
	else:
            azOffset = intDict[jdateid]['azInterp']
            elOffset = intDict[jdateid]['elInterp']
        ofile.write ("%d,%6.3f,%6.3f\n" % (obsnum, float(azOffset),float(elOffset)))
    
    ofile.close()


class LMTPointingTextFile:
    def __init__(self,filename):
        from csv import reader
        obsnums = []
        azOffset = []
        elOffset = []
        with open(filename, 'r') as csvfile:
            pReader = reader(csvfile, delimiter = ',', quotechar = "#")
            for irow in pReader:
                print irow
                obsnums.append(long(irow[0]))
                azOffset.append(float(irow[1]))
                elOffset.append(float(irow[2]))
        self.az = dict(zip(obsnums, azOffset))
        self.el = dict(zip(obsnums, elOffset))        
    
    
    
    
    
