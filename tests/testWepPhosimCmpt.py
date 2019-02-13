import os
import shutil
import numpy as np
import unittest

from lsst.ts.wep.Utility import FilterType

from lsst.ts.phosim.SkySim import SkySim

from lsst.ts.aoclcSim.Utility import getModulePath
from lsst.ts.aoclcSim.WepPhosimCmpt import WepPhosimCmpt


class TestWepPhosimCmpt(unittest.TestCase):
    """ Test the WepPhosimCmpt class."""

    def setUp(self):

        phosimDir = "phosimDir"
        self.wepPhosimCmpt = WepPhosimCmpt(phosimDir)

        # Set the output directories
        self.outputDir = os.path.join(getModulePath(), "tests", "tmp")
        self.outputImgDir = os.path.join(self.outputDir, "img")

        self.wepPhosimCmpt.setOutputDir(self.outputDir)
        self.wepPhosimCmpt.setOutputImgDir(self.outputImgDir)

        # Set the telescope survey parameters 
        obsId = 9006000
        filterType = FilterType.REF
        boresight = (0.2, 0.3)
        zAngleInDeg = 27.0912
        rotAngInDeg = np.rad2deg(-1.2323)
        mjd = 59552.3
        self.wepPhosimCmpt.setSurveyParam(obsId=obsId, filterType=filterType,
                                          boresight=boresight,
                                          zAngleInDeg=zAngleInDeg,
                                          rotAngInDeg=rotAngInDeg, mjd=mjd)

    def tearDown(self):

        shutil.rmtree(self.outputDir)

    def testGetOutputDir(self):

        outputDir = self.wepPhosimCmpt.getOutputDir()
        self.assertEqual(outputDir, self.outputDir)

    def testSetOutputDir(self):

        outputDir = os.path.join(self.outputDir, "testOutputDir")
        self.wepPhosimCmpt.setOutputDir(outputDir)

        self.assertTrue(self._isDirExists(outputDir))
        self.assertEqual(self.wepPhosimCmpt.getOutputDir(), outputDir)

    def _isDirExists(self, dirPath):

        return os.path.exists(dirPath)

    def testGetOutputImgDir(self):

        outputImgDir = self.wepPhosimCmpt.getOutputImgDir()
        self.assertEqual(outputImgDir, self.outputImgDir)

    def testSetOutputImgDir(self):

        outputImgDir = os.path.join(self.outputDir, "testOutputImgDir")
        self.wepPhosimCmpt.setOutputImgDir(outputImgDir)

        self.assertTrue(self._isDirExists(outputImgDir))
        self.assertEqual(self.wepPhosimCmpt.getOutputImgDir(), outputImgDir)

    def testGetSeedNum(self):

        self.assertEqual(self.wepPhosimCmpt.getSeedNum(), 0)

    def testSetSeedNum(self):

        seedNum = 3
        self.wepPhosimCmpt.setSeedNum(seedNum)

        self.assertEqual(self.wepPhosimCmpt.getSeedNum(), seedNum)

    def testGetPhosimParam(self):

        phosimParam = self.wepPhosimCmpt.getPhosimParam()

        self.assertEqual(phosimParam["numPro"], 1)
        self.assertEqual(phosimParam["e2ADC"], 1)

    def testSetPhosimParam(self):

        numPro = 3
        e2ADC = 0
        self._setPhosimParam(numPro, e2ADC)

        phosimParam = self.wepPhosimCmpt.getPhosimParam()

        self.assertEqual(phosimParam["numPro"], numPro)
        self.assertEqual(phosimParam["e2ADC"], e2ADC)

    def testSetPhosimParamWithWrongValue(self):

        numPro = -1
        e2ADC = 2
        self._setPhosimParam(numPro, e2ADC)

        phosimParam = self.wepPhosimCmpt.getPhosimParam()

        self.assertEqual(phosimParam["numPro"], 1)
        self.assertEqual(phosimParam["e2ADC"], 1)

    def _setPhosimParam(self, numPro, e2ADC):

        self.wepPhosimCmpt.setPhosimParam(numPro=numPro, e2ADC=e2ADC)

    def testGetOpdArgsAndFilesForPhoSim(self):

        self.wepPhosimCmpt.addOpdFieldXYbyDeg(0, 0)

        instFileName="opd.inst"
        argString = self.wepPhosimCmpt.getOpdArgsAndFilesForPhoSim(
                                                    instFileName=instFileName)

        self.assertEqual(len(argString), 240)
        self.assertEqual(self._getNumOfFileInFolder(self.outputDir), 11)

        instFilePath = os.path.join(self.outputDir, instFileName)
        numOfLine = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLine, 59)

    def testGetComCamOpdArgsAndFilesForPhoSim(self):

        instFileName="opd.inst"
        argString = self.wepPhosimCmpt.getComCamOpdArgsAndFilesForPhoSim(
                                                    instFileName=instFileName)

        self.assertEqual(len(argString), 240)
        self.assertEqual(self._getNumOfFileInFolder(self.outputDir), 11)

        instFilePath = os.path.join(self.outputDir, instFileName)
        numOfLine = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLine, 67)

    def _getNumOfFileInFolder(self, folder):
        
        return len([name for name in os.listdir(folder) 
                   if os.path.isfile(os.path.join(folder, name))])

    def _getNumOfLineInFile(self, filePath):

        with open(filePath, "r") as file:
            return sum(1 for line in file.readlines())

    def testGetStarArgsAndFilesForPhoSim(self):

        skySim = SkySim()
        skySim.addStarByRaDecInDeg(0, 0.1, 0.2, 5.0)

        instFileName="star.inst"
        argString = self.wepPhosimCmpt.getStarArgsAndFilesForPhoSim(
                                            skySim, instFileName=instFileName)

        self.assertEqual(len(argString), 243)
        self.assertEqual(self._getNumOfFileInFolder(self.outputDir), 11)

        instFilePath = os.path.join(self.outputDir, instFileName)
        numOfLine = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLine, 63)

    def testCalcOpdPssn(self):

        pssnList = self._getPssnListOfComCam()

        ansData = self._getAnsDataOfComCam()
        ansPssn = ansData[0, 0:9]

        pssnArray = np.array(pssnList)
        delta = np.sum(np.abs(pssnArray - ansPssn))
        self.assertLess(delta, 1e-10)

    def _getPssnListOfComCam(self):

        self._setOutputImgDir()
        pssnList = self.wepPhosimCmpt.calcComCamOpdPssn()[0]

        return pssnList

    def _setOutputImgDir(self):

        opdFileDir = self._getOpdFileDirOfComCam()
        self.wepPhosimCmpt.setOutputImgDir(opdFileDir)

    def _getOpdFileDirOfComCam(self):

        opdFileDir = os.path.join(getModulePath(), "tests", "testData",
                                  "comcamOpdFile", "iter0")

        return opdFileDir

    def _getAnsDataOfComCam(self):

        opdFileDir = self._getOpdFileDirOfComCam()
        ansFilePath = os.path.join(opdFileDir, "sim7_iter0_PSSN.txt")
        ansData = np.loadtxt(ansFilePath)

        return ansData

    def testCalcComCamGQeffFwhm(self):

        pssnList = self._getPssnListOfComCam()
        gqEffFwhm = self.wepPhosimCmpt.calcComCamOpdEffFwhm(pssnList)[1]

        ansData = self._getAnsDataOfComCam()
        ansGqEffFwhm = ansData[1, -1]

        delta = np.abs(ansGqEffFwhm - gqEffFwhm)
        self.assertLess(delta, 1e-10)

    def testMapOpdToZk(self):

        # Calculate the Zk from the OPD 
        self._setOutputImgDir()
        opdZkData = self.wepPhosimCmpt.mapOpdToZk()

        # Get the answer data
        opdFileDir = self._getOpdFileDirOfComCam()
        ansDataFilePath = os.path.join(opdFileDir, "sim7_iter0_opd.zer")
        ansOpdZkData = np.loadtxt(ansDataFilePath)[:, 3:]

        # Do the assertion
        delta = np.sum(np.abs(opdZkData - ansOpdZkData))
        self.assertLess(delta, 1e-10)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
