from aztec import *

def gaussian(x, c0,c1,x0,s0):
        return c0+ c1*npy.exp(-0.5*((x-x0)/s0)**2.)

def gaussian2d ((x,y), c0, c1, x0,y0, sx, sy, theta):
    x = npy.array(x)
    y = npy.array(y)
    xp = (x-x0)*npy.cos(theta) - (y-y0)*npy.sin(theta)
    yp = (x-x0)*npy.sin(theta) + (y-y0)*npy.cos(theta)
    expn = (xp/sx)**2.0 + (yp/sy)**2.0
    ff = c0 + c1*npy.exp(-0.5*expn)
    return ff

def gaussian2dplane ((x,y), c0, c1, x0,y0, sx, sy, theta,ax,ay):
    x = npy.array(x)
    y = npy.array(y)
    xp = (x-x0)*npy.cos(theta) - (y-y0)*npy.sin(theta)
    yp = (x-x0)*npy.sin(theta) + (y-y0)*npy.cos(theta)
    expn = (xp/sx)**2.0 + (yp/sy)**2.0
    ff = c0 + c1*npy.exp(-0.5*expn) +ax*xp + ay*yp
    return ff

def gaussian2dquad ((x,y), c0, c1, x0,y0, sx, sy, theta,ax,ay,bxy,bx2,by2):
    x = npy.array(x)
    y = npy.array(y)
    xp = (x-x0)*npy.cos(theta) - (y-y0)*npy.sin(theta)
    yp = (x-x0)*npy.sin(theta) + (y-y0)*npy.cos(theta)
    expn = (xp/sx)**2.0 + (yp/sy)**2.0
    ff = c0 + c1*npy.exp(-0.5*expn) +ax*xp + ay*yp + bxy*xp*yp + bx2*xp**2 + by2*yp**2
    return ff


def fitgaussian(data,x):
    from scipy.optimize import curve_fit
    pars,cov =curve_fit(gaussian, data, x)
    return  pars


def fitgaussian2d (data, x,y, sigma = None, iguess=None):

    from scipy.optimize import curve_fit
    
    xx, yy = npy.meshgrid(x,y)
    xx = xx.transpose()
    yy = yy.transpose()
    if iguess is None:
        iguess = [0.0, npy.max(data), 0.0, 0.0, 6.0/3600.0, 6.0/3600.0, 0.0]
    try:
        if not sigma is None:
            popt, pcov = curve_fit(gaussian2d,(xx.ravel(),yy.ravel()), data.ravel(), p0 = iguess, sigma=sigma.ravel(), maxfev=100000)
        else:
            popt, pcov = curve_fit(gaussian2d,(xx.ravel(),yy.ravel()), data.ravel(), p0 = iguess, maxfev=100000)
        return True,popt
    except Exception, err:
        print "Warning, this was a really bad fit: ", str(err)
        return False,-1.0*npy.ones(len(iguess))


def gaussian2dres (parameters, x,y, data, sigma):
    
    model = gaussian2d((x,y), parameters['dc'].value, parameters['flux'].value, parameters['xc'].value,
                       parameters['yc'].value, parameters['bx'].value, parameters['by'].value,
                       parameters['theta'].value)
    if sigma is None:
        sigma = npy.ones(model.shape)
    return (data-model)/sigma

