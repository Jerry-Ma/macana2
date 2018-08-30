from aztec import *
import os
from os import path as path
from fnmatch import fnmatch
from scipy.io.netcdf import netcdf_file as NetCDFFile

def file_search (pattern, directory=None, recursive = False):
    if directory is None:
        directory = "./"
    
    result = []
    
    if not path.exists(directory):
        raise Exception("aztec.util.file_search. Not such directory")
    
    for ifile in os.listdir(directory):
        if fnmatch(ifile,pattern):
            result.append(ifile)
    
    result.sort()
    return result

def findDataBounds (data):
    """ Find the bounding box [x1:x2,y1:y2] where data is non-zero
    """
    x1 = 0
    y1 = 0
    x2 = data.shape[0]-1
    y2 = data.shape[1]-1
    
    while npy.sum(data[x1,:]) == 0:
        x1+=1
    while npy.sum(data[x2,:]) == 0:
        x2-=1
    while npy.sum(data[:,y1]) == 0:
        y1+=1
    while npy.sum(data[:,y2]) == 0:
        y2-=1

    res ={'data': data[x1:x2,y1:y2].copy(), 'xi':x1, 'xe':x2, 'yi':y1, 'ye':y2}
    return res

def getJulDatefromNc (ncfile):
    #from netCDF4 import Dataset as NetCDFFile
    
    from astropy.time import Time
    from datetime import datetime
    nc = NetCDFFile (ncfile)
    
    retval = -1
    try:
        utdate = getattr(nc, 'date')
        uttime = getattr(nc, 'start_time')
        utdate = utdate.replace("\"", "")
        utString = "%s %s" %(utdate, uttime)
        t = Time(datetime.strptime(utString, "%m/%d/%Y %H:%M:%S"), scale="utc", format = "datetime")
        retval =  t.jd
    except Exception:
        #try:
         utdate = nc.variables['Header.TimePlace.UTDate'].getValue()
         t = Time (float(utdate), scale = "utc", format = "decimalyear")
         retval =  t.jd
        #except Exception:
            
    
    nc.close()
    return retval



def getAztecSeason (juldate):
    from datetime import datetime   
    from astropy.time import Time 
    
    dates = [
             [Time(datetime( 2005, 5, 1, 0, 0, 0), scale="utc", format = "datetime").jd,'05A'],
             [Time(datetime( 2005,10, 1, 0, 0, 0), scale="utc", format = "datetime").jd,'05B'],
             [Time(datetime( 2007,4, 1,  0, 0, 0), scale="utc", format = "datetime").jd,'ASTE07'],
             [Time(datetime( 2008,6, 1,  0, 0, 0), scale="utc", format = "datetime").jd,'ASTE08'],
             [Time(datetime( 2011,1, 1,  0, 0, 0), scale="utc", format = "datetime").jd,'LMT11A'],
             [Time(datetime( 2013,1, 1,  0, 0, 0), scale="utc", format = "datetime").jd,'LMT13A'],
             [Time(datetime( 2014,9, 1,  0, 0, 0), scale="utc", format = "datetime").jd,'LMTES3'],
             [Time(datetime( 2015,9, 1,  0, 0, 0), scale="utc", format = "datetime").jd,'LMTES4']
    ]

    season = "UNKNOWN"
    for i in dates[-1:0:-1]:
        if juldate > i[0]:
            season = i[1]
            break
    return season
        
def getObsNum(ncfile):


    nc = NetCDFFile (ncfile)
    try:    
        onum = nc.variables['Header.Dcs.ObsNum'].getValue()
    except:
        raise Exception("This does not seems to be a LMT file. Use ASTE functions instead")
    
    nc.close()
    return onum


def getOffsetFromMap(ncfile):


    nc = NetCDFFile (ncfile)
    try:  
        signalv = nc.variables['signal']  
        azOffset = signalv.offset_x
        elOffset = signalv.offset_y
    except KeyError:
        raise Exception("This does not seems to be a LMT file. Use ASTE functions instead")
    
    nc.close()
    return azOffset,elOffset

def isvalidObs(ncfile):


    nc = NetCDFFile (ncfile)
    try:  
        valid = nc.variables['Header.ScanFile.Valid'].getValue()
        if valid >0:
            return True
    except KeyError:
        print "Valid flag not found in nc file %s" % ncfile
        pass
    
    return False
