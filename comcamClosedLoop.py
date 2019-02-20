import os
import numpy as np

from lsst.ts.wep.Utility import FilterType, runProgram
from lsst.ts.wep.ctrlIntf.RawExpData import RawExpData
from lsst.ts.phosim.SkySim import SkySim

from lsst.ts.aoclcSim.Utility import getModulePath, getPhoSimPath
from lsst.ts.aoclcSim.WepPhosimCmpt import WepPhosimCmpt
from lsst.ts.aoclcSim.OfcCmpt import OfcCmpt


def _getComCamSensorNameList():

    sensorNameList = ["R22_S00", "R22_S01", "R22_S02", "R22_S10", "R22_S11",
                      "R22_S12", "R22_S20", "R22_S21", "R22_S22"]

    return sensorNameList


def _makeCalibs(outputDir, sensorNameList):

    fakeFlatDirName = "fake_flats"
    fakeFlatDir = os.path.join(outputDir, fakeFlatDirName)
    _makeDir(fakeFlatDir)

    detector = " ".join(sensorNameList)
    _genFakeFlat(fakeFlatDir, detector)


def _makeDir(directory):

    if (not os.path.exists(directory)):
        os.makedirs(directory)


def _genFakeFlat(fakeFlatDir, detector):
    
    currWorkDir = os.getcwd()

    os.chdir(fakeFlatDir)
    _makeFakeFlat(detector)
    os.chdir(currWorkDir)


def _makeFakeFlat(detector):

    command = "makeGainImages.py"
    argstring = "--detector_list %s" % detector
    runProgram(command, argstring=argstring)


def _prepareWepPhosimCmpt(phosimDirPath, filterType, rotAngInDeg, numPro):

    wepPhosimCmpt = WepPhosimCmpt(phosimDirPath)

    # Set the telescope survey parameters 
    boresight = (0, 0)
    zAngleInDeg = 27.0912
    mjd = 59552.3
    wepPhosimCmpt.setSurveyParam(filterType=filterType, boresight=boresight,
                                 zAngleInDeg=zAngleInDeg,
                                 rotAngInDeg=rotAngInDeg, mjd=mjd)

    # Set the PhoSim parameters
    wepPhosimCmpt.setPhosimParam(numPro=numPro, e2ADC=1)

    # Set the seed number for M1M3 surface
    seedNum = 6
    wepPhosimCmpt.setSeedNum(seedNum)

    return wepPhosimCmpt


def _prepareOfcCmpt(filterType, rotAngInDeg):

    ofcCmpt = OfcCmpt()
    ofcCmpt.setFilter(filterType)
    ofcCmpt.setRotAng(rotAngInDeg)

    return ofcCmpt


def _prepareSkySim(opdMetr):

    skySim = SkySim()

    starId = 0
    mag = 15
    raInDegList = opdMetr.fieldX
    declInDegList = opdMetr.fieldY
    for raInDeg, declInDeg in zip(raInDegList, declInDegList):
        skySim.addStarByRaDecInDeg(starId, raInDeg, declInDeg, mag)
        starId += 1

    return skySim


def main(phosimDirPath, iterNum, numPro):

    # Prepate the calibration products
    baseOutputDir = os.path.join(getModulePath(), "output")
    sensorNameList = _getComCamSensorNameList()
    _makeCalibs(baseOutputDir, sensorNameList)

    # Survey parameters
    filterType = FilterType.REF
    rotAngInDeg = 0.0

    # Prepare the components
    wepPhosimCmpt = _prepareWepPhosimCmpt(phosimDirPath, filterType,
                                          rotAngInDeg, numPro)
    ofcCmpt = _prepareOfcCmpt(filterType, rotAngInDeg)

    # Set the telescope state to be the same as the OFC
    state0 = ofcCmpt.getState0()
    wepPhosimCmpt.setDofInUm(state0)

    # Do the iteration
    obsId = 9006000
    opdZkFileName = "opd.zer"
    opdPssnFileName = "PSSN.txt"
    outputDirName = "pert"
    outputImgDirName = "img"
    iterDefaultDirName = "iter"
    dofInUmFileName = "dofPertInNextIter.mat"
    skyInfoFileName = "skyComCamInfo.txt"
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

        # Get the GQ effective FWHM from file
        gqEffFwhm = wepPhosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
        print(gqEffFwhm)

        # Prepare the faked sky according to the OPD field positions
        skySim = _prepareSkySim(wepPhosimCmpt.getOpdMetr())

        # Output the sky information. This is for the WepCmpt to use.
        outputSkyInfoFilePath = os.path.join(outputDir, skyInfoFileName)
        skySim.exportSkyToFile(outputSkyInfoFilePath)

        # Assign the entra- and intra-focal observation Id
        extraObsId = obsId + 1
        intraObsId = obsId + 2

        # Generate the defocal images
        argStringList = wepPhosimCmpt.getComCamStarArgsAndFilesForPhoSim(
            extraObsId, intraObsId, skySim, simSeed=1000,
            cmdSettingFileName="starDefault.cmd",
            instSettingFileName="starSingleExp.inst")
        for argString in argStringList:
            wepPhosimCmpt.runPhoSim(argString)

        # Repackage the images
        wepPhosimCmpt.repackageComCamImgFromPhoSim()

        # Get the PSSN from file
        pssn = wepPhosimCmpt.getOpdPssnFromFile(opdPssnFileName)
        print(pssn)

        # Set the gain value in OfcCmpt by pssn
        ofcCmpt.setGainByPSSN(pssn, sensorNameList)

        # Get the OPD zk from file
        opdZkData = wepPhosimCmpt.getZkFromFile(opdZkFileName)

        # Calculate the new DOF by OFC component
        dofInUm = ofcCmpt.calcAggDofForPhoSim(opdZkData, sensorNameList)

        # Set the new DOF to wepPhosimCmpt
        wepPhosimCmpt.setDofInUm(dofInUm)

        # Save the DOF file
        wepPhosimCmpt.saveDofInUmFileForNextIter(
                                    dofInUm, dofInUmFileName=dofInUmFileName)

        # Add the observation ID by 10 for the next iteration
        obsId += 10


if __name__ == "__main__":

    # PhoSim directory
    phosimDirPath = getPhoSimPath()

    # Processor Number
    numPro = 1

    # Iteration number
    iterNum = 1

    main(phosimDirPath, iterNum, numPro)
