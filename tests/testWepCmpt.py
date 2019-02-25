import os
import shutil
import numpy as np
import unittest

from lsst.ts.wep.Utility import FilterType, runProgram
from lsst.ts.wep.WepController import WepController
from lsst.ts.wep.ctrlIntf.RawExpData import RawExpData

from lsst.ts.aoclcSim.Utility import getModulePath
from lsst.ts.aoclcSim.WepCmpt import WepCmpt


class TestWepCmpt(unittest.TestCase):
    """ Test the WepCmpt class."""

    def setUp(self):

        self.outputDir = os.path.join(getModulePath(), "tests", "tmp")
        self._makeDir(self.outputDir)

        isrDirName = "input"
        isrDir = os.path.join(self.outputDir, isrDirName)
        self._makeDir(isrDir)

        self.wepCmpt = WepCmpt(isrDir)

        # Set the survey paramters
        self.wepCmpt.setFilter(FilterType.REF)
        self.wepCmpt.setBoresight(0.0, 0.0)
        self.wepCmpt.setRotAng(0.0)

    def _makeDir(self, newDir):

        os.makedirs(newDir, exist_ok=True)

    def tearDown(self):

        self.wepCmpt.disconnect()
        shutil.rmtree(self.outputDir)

    def testGetWepController(self):

        wepCntlr = self.wepCmpt.getWepController()
        self.assertTrue(isinstance(wepCntlr, WepController))

    def testGetFilter(self):

        filterType = self.wepCmpt.getFilter()
        self.assertEqual(filterType, FilterType.REF)

    def testSetFilter(self):

        filterType = FilterType.R
        self.wepCmpt.setFilter(filterType)

        self.assertEqual(self.wepCmpt.getFilter(), filterType)

    def testGetBoresight(self):

        raInDeg, decInDeg = self.wepCmpt.getBoresight()
        self.assertEqual(raInDeg, 0.0)
        self.assertEqual(decInDeg, 0.0)

    def testSetBoresight(self):

        raInDeg = 10.0
        decInDeg = 20.0
        self.wepCmpt.setBoresight(raInDeg, decInDeg)

        raInDegInWepCmpt, decInDegInWepCmpt = self.wepCmpt.getBoresight()
        self.assertEqual(raInDegInWepCmpt, raInDeg)
        self.assertEqual(decInDegInWepCmpt, decInDeg)

    def testGetRotAng(self):

        rotAngInDeg = self.wepCmpt.getRotAng()
        self.assertEqual(rotAngInDeg, 0.0)

    def testSetRotAng(self):

        rotAngInDeg = 10.0
        self.wepCmpt.setRotAng(rotAngInDeg)

        self.assertEqual(self.wepCmpt.getRotAng(), rotAngInDeg)

    def testIngestCalibs(self):

        sensorNameList = ["R22_S11"]
        fakeFlatDir = self._makeCalibs(self.outputDir, sensorNameList)

        numOfFile = self._getNumOfFileInFolder(fakeFlatDir)
        self.assertEqual(numOfFile, 6)

        self.wepCmpt.ingestCalibs(fakeFlatDir)

        numOfFile = self._getNumOfFileInFolder(fakeFlatDir)
        self.assertEqual(numOfFile, 0)

    def _makeCalibs(self, outputDir, sensorNameList):

        fakeFlatDirName = "fake_flats"
        fakeFlatDir = os.path.join(self.outputDir, fakeFlatDirName)
        self._makeDir(fakeFlatDir)

        detector = " ".join(sensorNameList)
        self._genFakeFlat(fakeFlatDir, detector)

        return fakeFlatDir

    def _genFakeFlat(self, fakeFlatDir, detector):

        currWorkDir = os.getcwd()

        os.chdir(fakeFlatDir)
        self._makeFakeFlat(detector)
        os.chdir(currWorkDir)

    def _makeFakeFlat(self, detector):

        command = "makeGainImages.py"
        argstring = "--detector_list %s" % detector
        runProgram(command, argstring=argstring)

    def _getNumOfFileInFolder(self, folder):

        return len([name for name in os.listdir(folder) 
                   if os.path.isfile(os.path.join(folder, name))])

    def testGetSkyFile(self):

        skyFile = self.wepCmpt.getSkyFile()
        self.assertEqual(skyFile, "")

    def testSetSkyFile(self):

        skyFile = "testSetSkyFile"
        self.wepCmpt.setSkyFile(skyFile)

        self.assertEqual(self.wepCmpt.getSkyFile(), skyFile)

    def testCalculateWavefrontErrorsComCam(self):

        # Make the calibration products and do the ingestion
        sensorNameList = ["R22_S11", "R22_S12"]
        fakeFlatDir = self._makeCalibs(self.outputDir, sensorNameList)
        self.wepCmpt.ingestCalibs(fakeFlatDir)

        # Set the skyFile
        repackagedDir = os.path.join(getModulePath(), "tests", "testData",
                                     "comcamRepackagedData")
        skyFilePath = os.path.join(repackagedDir, "skyComCamInfo.txt")
        self.wepCmpt.setSkyFile(skyFilePath)

        # Collect the wavefront data
        intraRawExpData = RawExpData()
        intraObsId = 9006002
        intraRawExpDir = os.path.join(repackagedDir, "intra")
        intraRawExpData.append(intraObsId, 0, intraRawExpDir)

        extraRawExpData = RawExpData()
        extraObsId = 9006001
        extraRawExpDir = os.path.join(repackagedDir, "extra")
        extraRawExpData.append(extraObsId, 0, extraRawExpDir)

        # Calculate the wavefront error
        wfErrMap = self.wepCmpt.calculateWavefrontErrorsComCam(intraRawExpData,
                                                               extraRawExpData)

        self.assertEqual(len(wfErrMap), 2)
        for wfErr in wfErrMap.values():
            self.assertEqual(wfErr.argmax(), 1)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
