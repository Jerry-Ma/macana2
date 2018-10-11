from aztec import *
from lxml.etree import Element, SubElement, tostring
from scipy.io.idl import readsav
from aztec.util import getAztecSeason, getJulDatefromNc
from scipy.io.netcdf import netcdf_file
from scipy.io.netcdf import netcdf_file as NetCDFFile
#from lxml.etree import tostring

def addVariable (topElement, name, value, units=None):
    
    var = SubElement (topElement, name)
    subEl = SubElement (var, "value")
    subEl.text = str(value)
    if not units is None:
        subEl2 = SubElement (var, "units")
        subEl2.text = str (units)

def boloInfo2Xml (filename, boloInfo):
    
    
    nbolo = len (boloInfo)
    top = Element('nBolos')
    xnbolo = SubElement(top,"value")
    xnbolo.text = str (nbolo)
    
    xmlFullString = tostring(top, pretty_print=True)
    
    for ibolo, i in zip(boloInfo,range(nbolo)):
        top = Element ("d%d"%i)
        addVariable (top, "name",ibolo['boloName'])
        addVariable (top, "ncdf_location", ibolo['id'])
        
        
        
        #Tau fit
        addVariable(top, "az_fwhm", ibolo['az_fwhm'], units = "arcseconds")
        addVariable(top, "el_fwhm", ibolo['el_fwhm'], units = "arcseconds")
        addVariable(top, "az_offset", ibolo['az_offset'], units = "arcseconds")
        addVariable(top, "el_offset", ibolo['el_offset'], units = "arcseconds")
        addVariable(top, "bologain", ibolo['bologain'], units = "mJy/nW")
        addVariable(top, "bologain_err", ibolo['bologain_err'], units = "mJy/nW")
        addVariable(top, "bolosens", ibolo['bolosens'], units = "mJy-rt(s)")
        addVariable(top, "offset_dc2tau", ibolo['tau_offset'])
        addVariable(top, "offset_err_dc2tau", ibolo['tau_offset_err'])
        addVariable(top, "slope_dc2tau", ibolo['tau_slope'])
        addVariable(top, "slope_err_dc2tau", ibolo['tau_slope_err'])
        addVariable(top, "quad_dc2tau", ibolo['tau_quad'])
        addVariable(top, "quad_err_dc2tau", ibolo['tau_quad_err'])
        addVariable(top, "offset_dc2responsivity", ibolo['res_offset'])
        addVariable(top, "offset_err_dc2responsivity", ibolo['res_offset_err'])
        addVariable(top, "slope_dc2responsivity", ibolo['res_slope'])
        addVariable(top, "slope_err_dc2responsivity", ibolo['res_slope_err'])
        addVariable (top, "goodflag", ibolo['valid'])
#         if ibolo['mean_el'] != -1:
#             addVariable(top, "mean_el", ibolo['mean_el'], units = "degrees")
        xmlFullString += tostring (top, pretty_print=True)       
    
    ff = open (filename, "w")
    ff.write(xmlFullString)
    ff.close()

def fcf_model (p, el, afs, ibolo):
    """
        This is a slightly modified versions of the idl fcf_model on IDL pipeline
    """
    el_model = p[0] + p[1]*el + p[2]*el**2. + p[3]*el**3. + p[4]*el**4. + p[5]*el**5.
    afs_model = p[6]*afs 
    return p[7+ibolo]*(el_model + afs_model)


def fix_nw_bolos (jdate, boloInfo):
    
    if jdate >2456966.05 and jdate < 2456967.12: #Bolometers in hex 5 were not working
        for ibolo in boloInfo:
            if ibolo['boloName'].find("h5") != -1:
                ibolo['valid']=0
                print "Removing bolometer %s. It was not working on date %d" %(ibolo['boloName'],jdate)

def interpolateParams (jdate, boloInfo, meanEl = None):
    season =getAztecSeason(jdate) 
    #Deal with non-uniform name tagging
    allBolo = False
    if season.find("ASTE") != -1:
        season2 = season.lower()
    elif season.find("LMT") != -1:
        season2 = "LMT"
        allBolo = True
    else:
        season2 = ""
        
    calpath = os.path.join(os.getenv("AZTEC_MACANA_PATH"),"calibration/")    
    beamparssave = os.path.join(calpath, "parameters_%s" %season, "aztec_%s_params_beammap.sav" % season2)
    pars = readsav(beamparssave)
    
    print "Using file %s for Julian Date %d" % (beamparssave, jdate)
    if "jd_all" in pars.keys():                     #Stuff for ASTE data
        ss = npy.argsort (pars['jd_all'])
        jd = pars['jd_all'][ss]
        azfwhm = pars['azfwhm'][ss,:].T
        azfwhmerr = pars['azfwhmerr'][ss,:].T
        elfwhm = pars['elfwhm'][ss,:].T
        elfwhmerr = pars['elfwhmerr'][ss,:].T
        azoff = pars['azoff'][ss,:].T
        azofferr = pars['azofferr'][ss,:].T
        eloff = pars['eloff'][ss,:].T
        elofferr = pars['elofferr'][ss,:].T
        gain = pars['gain'][ss,:].T
        gainerr = pars['gainerr'][ss,:].T
        sens = pars['sens'][ss,:].T
        sources = pars ['sources'][ss]
        p=None
        meanel=None
        shiftStart = pars['shift_jd_low']
        shiftEnd = pars ['shift_jd_high']
        ishift = npy.where (npy.logical_and(jdate >= shiftStart, jdate <=shiftEnd))
        shiftStart = shiftStart[ishift]
        shiftEnd = shiftEnd[ishift]
        #Now we get the mean values
        meanAzFwhm = pars ['mean_azfwhm']
        meanAzOff = pars [ 'mean_azoff']
        meanElFwhm = pars [ 'mean_elfwhm']
        meanElOff = pars [ 'mean_eloff']
        meanGain = pars [ 'mean_gain']
        meanGainErr = pars [ 'mean_gain_obs']
        meanSens = pars [ 'mean_sens']
        meanSensErr = pars [ 'mean_senssub']
 
    elif "jd_beammap" in pars.keys():
        ss = npy.argsort (pars['jd_beammap'])
        jd = pars['jd_beammap'][ss]
        azfwhm = pars['az_fwhm_array'][:,ss]
        azfwhmerr = pars['det_std_az_fwhm']
        elfwhm = pars['el_fwhm_array'][:,ss]
        elfwhmerr = pars['det_std_el_fwhm']
        azoff = pars['az_offset_array'][:,ss]
        azofferr = pars['det_std_az_off']
        eloff = pars['el_offset_array'][:,ss]
        elofferr = pars['det_std_el_off']
        gain = None
        gainerr = None
        sens = pars['sens_array'][:,ss]
        sources = pars ['sources'][ss]
        p= pars['p']
        meanel = pars['el_beammap'][ss]
        shiftStart = npy.floor(jdate) + 5./24.
        shiftEnd = shiftStart+1
        #Now we get the mean values
        meanAzFwhm = pars ['det_med_az_fwhm']
        meanAzOff = pars [ 'det_med_az_off']
        meanElFwhm = pars [ 'det_med_el_fwhm']
        meanElOff = pars [ 'det_med_el_off']
        meanSens = pars [ 'det_med_sens']
        meanAzFwhmErr = pars [ 'det_std_az_fwhm']
        meanSensErr = pars [ 'det_std_sens']
        
    else:
        raise Exception("No beam params information")
    
    
    #bracktec    
    posn = npy.abs(jdate-jd).argmin()
    if jdate < jd.min():
        bin0, bin1 = 0,0
    elif  jdate >= jd.max():
        bin0, bin1 = jd.shape[0]-1, jd.shape[0]-1
    elif jd[posn] == jdate:
        bin0, bin1 = posn,posn
    elif jd[posn] < jdate:
        bin0 = posn
        bin1 = posn+1
    else:
        bin0 = posn-1
        bin1 = posn
  

    #Now check if the selected beammaps were taken on the same or shift
    useMean= False
    
    if jd[bin0] < shiftStart and jd[bin1] > shiftEnd:
        useMean = True
    elif jd[bin0] < shiftStart and jd[bin1]<= shiftEnd:
        bin0 = bin1
    elif jd[bin0] >=shiftStart and jd[bin1]>shiftEnd:
        bin1 = bin0
    else:
        pass
    
#    import ipdb; ipdb.set_trace()
#     
#     if not npy.isfinite(azfwhmerr[bin0]) or not npy.isfinite(azfwhmerr[bin1]):
#         useMean = True
 
    
    fix_nw_bolos(jdate, boloInfo)
    
    nbolo = len(boloInfo)
    j = 0
    
    if not useMean:
        for ibolo,i in zip(boloInfo, range(nbolo)):
            if ibolo['valid'] == 1:
                if allBolo:     #LMT info has all bolometers (working and not working)
                    k=i
                else:           #ASTE has only the working bolometer info
                    k=j
                ibolo['az_fwhm'] = npy.interp(jdate, [jd[bin0],jd[bin1]], [azfwhm[k,bin0],azfwhm[k,bin1]])
                ibolo['el_fwhm'] = npy.interp(jdate, [jd[bin0],jd[bin1]], [elfwhm[k,bin0],elfwhm[k,bin1]])
                ibolo['az_offset'] = npy.interp(jdate, [jd[bin0],jd[bin1]], [azoff[k,bin0],azoff[k,bin1]])
                ibolo['el_offset'] = npy.interp(jdate, [jd[bin0],jd[bin1]], [eloff[k,bin0],eloff[k,bin1]])
                ibolo['bolosens'] = npy.interp(jdate, [jd[bin0],jd[bin1]], [sens[k,bin0],sens[k,bin1]])
                if gain is None:
                    ibolo['bologain'] = fcf_model (p,meanEl,azfwhmerr[k],j)
                    ibolo['bologain_err'] = npy.abs(ibolo['bologain'] * 0.1)     #Temporary hack as in IDL pipeline
                    ibolo ['mean_el'] = meanEl
                else:
                    ibolo['bologain'] = npy.interp(jdate, [jd[bin0],jd[bin1]], [gain[k,bin0],gain[k,bin1]])
                    ibolo['bologain_err'] = npy.interp(jdate, [jd[bin0],jd[bin1]], [gainerr[k,bin0],gainerr[k,bin1]])
                    ibolo['mean_el'] = -1
                j+=1
            else:
                ibolo['az_fwhm'] = npy.nan
                ibolo['el_fwhm'] = npy.nan
                ibolo['az_offset'] = npy.nan
                ibolo['el_offset'] = npy.nan
                ibolo['bolosens'] = npy.nan
                ibolo['bologain'] = npy.nan
                ibolo['bologain_err'] = npy.nan
                if gain is None:
                    ibolo['mean_el'] = npy.nan
                else:
                    ibolo['mean_el'] = -1
    else:
        print "Use mean bolometer parameters"
        for ibolo,i in zip(boloInfo, range(nbolo)):
            if ibolo['valid'] == 1:
                if allBolo:     #LMT info has all bolometers (working and not working)
                    k=i
                else:           #ASTE has only the working bolometer info
                    k=j
                ibolo['az_fwhm'] = meanAzFwhm[k]
                ibolo['el_fwhm'] = meanElFwhm[k]
                ibolo['az_offset'] = meanAzOff[k]
                ibolo['el_offset'] = meanElOff[k]
                ibolo['bolosens'] = meanSens[k]
                if gain is None:
                    map_afs = meanAzFwhmErr[k]
                    ibolo['bologain'] = fcf_model (p,meanEl,map_afs,j)
                    ibolo['bologain_err'] = ibolo['bologain'] * 0.1     #Temporary hack as in IDL pipeline
                    ibolo ['mean_el'] = meanEl
                else:
                    ibolo['bologain'] = meanGain[k]
                    ibolo['bologain_err'] = meanGainErr[k]
                    ibolo['mean_el'] = -1
                j+=1
            else:
                ibolo['az_fwhm'] = npy.nan
                ibolo['el_fwhm'] = npy.nan
                ibolo['az_offset'] = npy.nan
                ibolo['el_offset'] = npy.nan
                ibolo['bolosens'] = npy.nan
                ibolo['bologain'] = npy.nan
                ibolo['bologain_err'] = npy.nan
                if gain is None:
                    ibolo['mean_el'] = npy.nan
                else:
                    ibolo['mean_el'] = -1

def createBstatsFile (ncfile, outFile = None):
    
    from scipy.interpolate import interp1d
    from xml.etree.ElementTree import Element, SubElement
    
    if outFile is None:
        outFile = ncfile.replace (".nc", ".bstats")
        
    jdate = getJulDatefromNc(ncfile)
    season =getAztecSeason(jdate) 
    
    if season.find("LMT")>-1:
        
        nc = NetCDFFile(ncfile)
        elSignal = npy.mean(nc.variables['Data.AztecBackend.TelElAct'][:])
        elSignal = elSignal.flatten()
        meanEl = npy.degrees(npy.median (elSignal))
        nc.close()
        namePrefix = "Data.AztecBackend."
    else:
        meanEl = None
        namePrefix = ""
    
    calpath = os.path.join(os.getenv("AZTEC_MACANA_PATH"),"calibration/")
    bololistfile = os.path.join(calpath, "parameters_%s" %season, "bololist.csv")
    tausave = os.path.join(calpath, "parameters_%s" %season, "fit_parameters_bolodc2tau_%s.sav" % season)    
    ressave = os.path.join(calpath, "parameters_%s" %season, "fit_parameters_bolodc2responsivity_%s.sav" % season)
    

    bololist = npy.loadtxt (bololistfile, delimiter = ",", comments = "#", dtype={'names': ('boloname', 'flags',),
                                                                                 'formats': ('|S10', 'i4')})
    
    nbolo = len (bololist)
    boloInfo = []
    jbolo = 1
    for ibolo in bololist:
        boloInfo.append({'boloName':namePrefix+ibolo[0], 'valid':ibolo[1], 'id':jbolo})
        jbolo+=1
    
    tauInfo = readsav(tausave)
    resInfo = readsav(ressave)
    
    for ibolo, i in zip(boloInfo, range(nbolo)):
        
        if "offset" in tauInfo.keys():
            ibolo['tau_offset'] = tauInfo['offset'][i]
            ibolo['tau_offset_err'] = tauInfo['offset_err'][i]
            ibolo['tau_slope'] = tauInfo['slope'][i]
            ibolo['tau_slope_err'] = tauInfo['slope_err'][i]
            ibolo['tau_quad'] = 0.0
            ibolo['tau_quad_err'] = 0.0
        elif "p0" in tauInfo.keys():
            ibolo['tau_offset'] = tauInfo['p0'][i]
            ibolo['tau_offset_err'] = tauInfo['p0_err'][i]
            ibolo['tau_slope'] = tauInfo['p1'][i]
            ibolo['tau_slope_err'] = tauInfo['p1_err'][i]
            ibolo['tau_quad'] = tauInfo['p2'][i]
            ibolo['tau_quad_err'] = tauInfo['p2_err'][i]         
        else:
            raise Exception("Cannot get Tau2dc information properly.")
        ibolo['res_offset'] = resInfo['offset'][i]
        ibolo['res_offset_err'] = resInfo['offset_err'][i]
        ibolo['res_slope'] = resInfo['slope'][i]
        ibolo['res_slope_err'] = resInfo['slope_err'][i]

#Now interpolate the gain values
    interpolateParams(jdate, boloInfo, meanEl=meanEl)
    boloInfo2Xml(outFile,boloInfo)     
    
