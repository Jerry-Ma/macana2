import os

def makeDirectoryStructure(mdir, pointing=False):
    macana_dirs = ["raw_data", "reduced_maps", "noise_maps", "coadded_maps"]
    for idir in macana_dirs:
        if not os.path.exists(os.path.join(mdir,idir)):
            os.mkdir(os.path.join(mdir,idir))
    if pointing:
        pointing_dirs = ["raw_data", "reduced_maps"]
        pdir = os.path.join(mdir, "pointing")
        if not os.path.exists(pdir):
            os.mkdir(pdir)
        for idir in macana_dirs:
            if not os.path.exists(os.path.join(pdir,idir)):
                os.mkdir(os.path.join(pdir,idir))
        

def runMacana(apfile, macanaPath = os.path.expanduser(os.environ['AZTEC_MACANA_PATH'])):
    import subprocess
    curDir = os.getcwd()
    apfile = os.path.abspath(apfile)
    command = macanaPath + "/macana"
    subprocess.check_call([command, apfile])
    os.chdir(curDir)