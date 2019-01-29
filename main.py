import os
import numpy as np

from lsst.ts.wep.Utility import FilterType

from lsst.ts.aoclcSim.Utility import getModulePath
from lsst.ts.aoclcSim.WepPhosimCmpt import WepPhosimCmpt


def _prepareWepPhosimCmpt(phosimDirPath):

    wepPhosimCmpt = WepPhosimCmpt(phosimDirPath)

    # Set the telescope survey parameters 
    filterType = FilterType.REF
    boresight = (0, 0)
    zAngleInDeg = 27.0912
    rotAngInDeg = 0.0
    mjd = 59552.3
    wepPhosimCmpt.setSurveyParam(filterType=filterType, boresight=boresight,
                                 zAngleInDeg=zAngleInDeg,
                                 rotAngInDeg=rotAngInDeg, mjd=mjd)

    # Set the PhoSim parameters
    wepPhosimCmpt.setPhosimParam(numPro=1, e2ADC=1)

    # Set the seed number for M1M3 surface
    seedNum = 6
    wepPhosimCmpt.setSeedNum(seedNum)

    return wepPhosimCmpt


def main(phosimDirPath, iterNum):

    # Prepare the components
    wepPhosimCmpt = _prepareWepPhosimCmpt(phosimDirPath)

    # Do the iteration
    obsId = 9006000
    baseOutputDir = os.path.join(getModulePath(), "output")
    for iterCount in range(iterNum):

        # Set the observation Id
        wepPhosimCmpt.setSurveyParam(obsId=obsId)

        # The iteration directory
        iterDirName = "iter%d" % iterCount

        # Set the output directory
        outputDir = os.path.join(baseOutputDir, iterDirName, "pert")
        wepPhosimCmpt.setOutputDir(outputDir)

        # Set the output image directory
        outputImgDir = os.path.join(baseOutputDir, iterDirName, "img")
        wepPhosimCmpt.setOutputImgDir(outputImgDir)

        # Update the DOF
        dofInUm = np.zeros(50)
        dofInUm[5] = wepPhosimCmpt.tele.getDefocalDisInMm() * 1e3
        if (iterCount > 0):
            wepPhosimCmpt.accDofInUm(dofInUm)

        # Generate the OPD image
        argString = wepPhosimCmpt.getComCamOpdArgsAndFilesForPhoSim()
        wepPhosimCmpt.runPhoSim(argString)

        # Calculate the PSSN
        pssnList = wepPhosimCmpt.calcOpdPssn()

        # Calculate the GQ effective FWHM
        gqEffFwhm = wepPhosimCmpt.calcComCamGQeffFwhm(pssnList)
        print(gqEffFwhm)

        # Add the observation ID by 1
        obsId += 1


if __name__ == "__main__":

    # PhoSim directory
    phosimDirPath = os.path.join(os.sep, "home", "lsst", "phosim_syseng4")

    # Iteration number
    iterNum = 2
    main(phosimDirPath, iterNum)
