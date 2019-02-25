import os
import numpy as np

from lsst.ts.wep.Utility import FilterType

from lsst.ts.aoclcSim.Utility import getModulePath, getPhoSimPath
from lsst.ts.aoclcSim.WepPhosimCmpt import WepPhosimCmpt
from lsst.ts.aoclcSim.OfcCmpt import OfcCmpt


def _prepareWepPhosimCmpt(phosimDirPath, filterType, rotAngInDeg):

    wepPhosimCmpt = WepPhosimCmpt(phosimDirPath)

    # Set the telescope survey parameters 
    boresight = (0, 0)
    zAngleInDeg = 27.0912
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


def _prepareOfcCmpt(filterType, rotAngInDeg):

    ofcCmpt = OfcCmpt()
    ofcCmpt.setFilter(filterType)
    ofcCmpt.setRotAng(rotAngInDeg)

    return ofcCmpt


def _getComCamSensorNameList():

    sensorNameList = ["R22_S00", "R22_S01", "R22_S02", "R22_S10", "R22_S11",
                      "R22_S12", "R22_S20", "R22_S21", "R22_S22"]
    return sensorNameList


def main(phosimDirPath, iterNum):

    # Survey parameters
    filterType = FilterType.REF
    rotAngInDeg = 0.0

    # Prepare the components
    wepPhosimCmpt = _prepareWepPhosimCmpt(phosimDirPath, filterType,
                                          rotAngInDeg)
    ofcCmpt = _prepareOfcCmpt(filterType, rotAngInDeg)

    # Set the telescope state to be the same as the OFC
    state0 = ofcCmpt.getState0()
    wepPhosimCmpt.setDofInUm(state0)

    # Get the sensor name of ComCam
    sensorNameList = _getComCamSensorNameList()

    # Do the iteration
    obsId = 9006000
    baseOutputDir = os.path.join(getModulePath(), "output")
    opdZkFileName = "opd.zer"
    opdPssnFileName = "PSSN.txt"
    outputDirName = "pert"
    outputImgDirName = "img"
    iterDefaultDirName = "iter"
    dofInUmFileName="dofPertInNextIter.mat"
    for iterCount in range(iterNum):

        # Set the observation Id
        wepPhosimCmpt.setSurveyParam(obsId=obsId)

        # The iteration directory
        iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

        # Set the output directory
        outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
        wepPhosimCmpt.setOutputDir(outputDir)

        # Set the output image directory
        outputImgDir = os.path.join(baseOutputDir, iterDirName,
                                    outputImgDirName)
        wepPhosimCmpt.setOutputImgDir(outputImgDir)

        # Generate the OPD image
        argString = wepPhosimCmpt.getComCamOpdArgsAndFilesForPhoSim()
        wepPhosimCmpt.runPhoSim(argString)

        # Analyze the OPD data
        wepPhosimCmpt.analyzeComCamOpdData(zkFileName=opdZkFileName,
                                           pssnFileName=opdPssnFileName)

        # Get the PSSN from file
        pssn = wepPhosimCmpt.getOpdPssnFromFile(opdPssnFileName)
        print(pssn)

        # Set the gain value in OfcCmpt by pssn
        ofcCmpt.setGainByPSSN(pssn, sensorNameList)

        # Get the GQ effective FWHM from file
        gqEffFwhm = wepPhosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
        print(gqEffFwhm)

        # Get the OPD zk from file
        opdZkData = wepPhosimCmpt.getZkFromFile(opdZkFileName)

        # Calculate the new DOF by OFC component
        dofInUm = ofcCmpt.calcAggDofForPhoSim(opdZkData, sensorNameList)

        # Set the new DOF to wepPhosimCmpt
        wepPhosimCmpt.setDofInUm(dofInUm)

        # Save the DOF file
        wepPhosimCmpt.saveDofInUmFileForNextIter(
                                    dofInUm, dofInUmFileName=dofInUmFileName)

        # Add the observation ID by 1
        obsId += 1


if __name__ == "__main__":

    # PhoSim directory
    phosimDirPath = getPhoSimPath()

    # Iteration number
    iterNum = 5

    main(phosimDirPath, iterNum)
