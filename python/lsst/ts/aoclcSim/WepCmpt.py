import os
import numpy as np

from lsst.ts.wep.Utility import getModulePath, FilterType, CamType, BscDbType, \
                                abbrevDectectorName

from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.WfEstimator import WfEstimator
from lsst.ts.wep.WepController import WepController

import lsst.ts.aoclcSim.Utility as aoclcUtil


class WepCmpt(object):

    STAR_RADIUS_IN_PIXEL = 63
    SPACING_COEFFICIENT = 2.5
    MAX_NEIGHBORING_STAR = 1

    # with 1.5 mm defocal distance
    DEFOCAL_DIS_IN_MM = 1.5
    DONUT_IMG_SIZE_IN_PIXEL = 160

    def __init__(self, isrDir):
        """Initialization of WEP component class.

        Parameters
        ----------
        isrDir : str
            ISR directory. This directory path will have the input and output
            that the data butler needs.
        """

        self.isrDir = isrDir
        self.raInDeg = 0.0
        self.decInDeg = 0.0
        self.rotSkyPos = 0.0
        self.skyFile = ""

        # Configure the WEP controller
        self.wepCntlr = self._configWepController(isrDir)

    def _configWepController(self, isrDir):
        """Configure the WEP controller.

        WEP: wavefront estimation pipeline.

        Parameters
        ----------
        isrDir : str
            ISR directory. This directory path will have the input and output
            that the data butler needs.

        Returns
        -------
        WepController
            Configured WEP controller.
        """

        dataCollector = CamDataCollector(isrDir)
        isrWrapper = CamIsrWrapper(isrDir)
        sourSelc = self._configSourceSelector()
        sourProc = self._configSourceProcessor()
        wfsEsti = self._configWfEstimator()
        wepCntlr = WepController(dataCollector, isrWrapper, sourSelc,
                                 sourProc, wfsEsti)

        return wepCntlr

    def _configSourceSelector(self):
        """Configue the source selector.

        Returns
        -------
        SourceSelector
            Configured source selector.
        """

        sourSelc = SourceSelector(CamType.ComCam, BscDbType.LocalDbForStarFile)

        # Set the criteria of neighboring stars
        sourSelc.configNbrCriteria(self.STAR_RADIUS_IN_PIXEL,
                                   self.SPACING_COEFFICIENT, 
                                   maxNeighboringStar=self.MAX_NEIGHBORING_STAR)

        # Connest the database
        dbAdress = os.path.join(aoclcUtil.getModulePath(), "tests", "testData",
                                "bsc.db3")
        sourSelc.connect(dbAdress)

        return sourSelc

    def _configSourceProcessor(self):
        """Configure the source processor.

        Returns
        -------
        SourceProcessor
            Configured source processor.
        """

        folderPath2FocalPlane = os.path.join(getModulePath(), "tests",
                                             "testData")
        sourProc = SourceProcessor(folderPath2FocalPlane=folderPath2FocalPlane)

        return sourProc

    def _configWfEstimator(self):
        """Configure the wavefront estimator.

        Returns
        -------
        WfEstimator
            Configured wavefront estimator.
        """

        instruFolderPath = os.path.join(self._getConfigDataPath(), "cwfs",
                                        "instruData")
        algoFolderPath = os.path.join(self._getConfigDataPath(), "cwfs",
                                      "algo")
        wfsEsti = WfEstimator(instruFolderPath, algoFolderPath)

        # Use the comcam to calculate the LSST central raft image
        wfsEsti.config(solver="exp", instName="comcam",
                       opticalModel="offAxis",
                       defocalDisInMm=self.DEFOCAL_DIS_IN_MM,
                       sizeInPix=self.DONUT_IMG_SIZE_IN_PIXEL, debugLevel=0)

        return wfsEsti

    def _getConfigDataPath(self):
        """Get the configuration data path.

        Returns
        -------
        str
            Configuration data path
        """

        configDataPath = os.path.join(getModulePath(), "configData")

        return configDataPath

    def getWepController(self):
        """Get the WEP controller.

        WEP: wavefront estimation pipeline.

        Returns
        -------
        WepController
            WEP controller.
        """

        return self.wepCntlr

    def disconnect(self):
        """Disconnect the database."""

        self.wepCntlr.sourSelc.disconnect()

    def setFilter(self, filterType):
        """Set the current filter.

        Parameters
        ----------
        filterType : FilterType
            The new filter configuration to use for WEP data processing.
        """

        self.wepCntlr.sourSelc.setFilter(filterType)

    def getFilter(self):
        """Get the current filter.

        Returns
        -------
        FilterType
            The current filter configuration to use for WEP data processing.
        """

        return self.wepCntlr.sourSelc.getFilter()

    def setBoresight(self, raInDeg, decInDeg):
        """Set the boresight (ra, dec) in degree from the pointing component.

        The cooridinate system of pointing component is the international
        cannabinoid research society (ICRS).

        Parameters
        ----------
        raInDeg : float
            Right ascension in degree. The value should be in (0, 360).
        decInDeg : float
            Declination in degree. The value should be in (-90, 90). 
        """

        self.raInDeg = raInDeg
        self.decInDeg = decInDeg

    def getBoresight(self):
        """Get the boresight (ra, dec) defined in the international
        cannabinoid research society (ICRS).

        Returns
        -------
        raInDeg : float
            Right ascension in degree. The value should be in (0, 360).
        decInDeg : float
            Declination in degree. The value should be in (-90, 90). 
        """

        return self.raInDeg, self.decInDeg

    def setRotAng(self, rotAngInDeg):
        """Set the camera rotation angle in degree from the camera rotator
        control system.

        Parameters
        ----------
        rotAngInDeg : float
            The camera rotation angle in degree (-90 to 90).
        """

        self.rotSkyPos = rotAngInDeg

    def getRotAng(self):
        """Get the camera rotation angle in degree defined in the camera
        rotator control system.

        Returns
        -------
        float
            The camera rotation angle in degree.
        """

        return self.rotSkyPos

    def ingestCalibs(self, calibsDir):
        """Ingest the calibration products.

        Parameters
        ----------
        calibsDir : str
            Calibration directory.
        """

        mapperFile = os.path.join(self.isrDir, "_mapper")
        if (not os.path.exists(mapperFile)):
            self.wepCntlr.dataCollector.genPhoSimMapper()

        calibFiles = os.path.join(calibsDir, "*")
        self.wepCntlr.dataCollector.ingestCalibs(calibFiles)

    def setSkyFile(self, skyFile):
        """Set the sky information file.

        This is a temporary function for WepCmpt to query the star file
        generated by WepPhosimCmpt. This function will be removed after we begin
        to integrate the bright star catalog database.

        Parameters
        ----------
        skyFile : str
            Sky information file.
        """

        self.skyFile = skyFile

    def getSkyFile(self):
        """Get the sky information file.

        Returns
        -------
        str
            Sky information file.
        """

        return self.skyFile

    def calculateWavefrontErrorsComCam(self, intraRawExpData, extraRawExpData):
        """Calculate the wavefront errors for ComCam.

        ComCam: commissioning camera.

        Parameters
        ----------
        intraRawExpData : RawExpData
            Intra-focal raw exposure data.
        extraRawExpData : RawExpData
            Extra-focal raw exposure data.

        Returns
        -------
        dict
            Calculated wavefront error. The dictionary key [str] is the
            abbreviated sensor name (e.g. R22_S11). The dictionary item
            [numpy.ndarray] is the averaged wavefront error (z4-z22) in nm. 
        """

        # Do the ingestion of images
        self._ingestImg(intraRawExpData)
        self._ingestImg(extraRawExpData)

        # Do the ISR
        fileName = "isr_config.py"
        isrFilePath = os.path.join(self.isrDir, fileName)
        self.wepCntlr.isrWrapper.config(doFlat=True, fileName=fileName)

        rerunName = "run1"
        self.wepCntlr.isrWrapper.doISR(self.isrDir, rerunName=rerunName)

        # Set Butler Inputs Path
        inputs = os.path.join(self.isrDir, "rerun", rerunName)
        self.wepCntlr.setPostIsrCcdInputs(inputs)

        # Get the target star by file
        self.wepCntlr.sourSelc.setObsMetaData(self.raInDeg, self.decInDeg,
                                              self.rotSkyPos)
        neighborStarMap, starMap, wavefrontSensors = \
            self.wepCntlr.sourSelc.getTargetStarByFile(self.skyFile, offset=0)

        # Get the sensor name list
        sensorNameList = list(wavefrontSensors)

        # Begin to calculate the wavefront error
        intraObsIdList = intraRawExpData.getVisit()
        extraObsIdList = extraRawExpData.getVisit()
        avgErrMapList = []
        for intraObsId, extraObsId in zip(intraObsIdList, extraObsIdList):

            # Get Post-ISR defocal image map
            obsIdList = [intraObsId, extraObsId]
            wfsImgMap = self.wepCntlr.getPostIsrImgMapByPistonDefocal(
                                                sensorNameList, obsIdList)

            # Get the donut map
            donutMap = self.wepCntlr.getDonutMap(
                        neighborStarMap, wfsImgMap, self.getFilter(),
                        doDeblending=False)

            # Calculate the wavefront error
            donutMap = self.wepCntlr.calcWfErr(donutMap)

            # Calculate the average wavefront error on single CCD
            avgErrMap = dict()
            for sensor, donutList in donutMap.items():
                avgErr = self.wepCntlr.calcAvgWfErrOnSglCcd(donutList)

                abbrevSensor = abbrevDectectorName(sensor)
                avgErrMap[abbrevSensor] = avgErr

            # Collect the average wavefront error map to the list
            avgErrMapList.append(avgErrMap)

        # Since there is only one pait at this moment, just return the first one.

        return avgErrMapList[0]

    def _ingestImg(self, rawExpData):
        """Ingest the images.

        Parameters
        ----------
        rawExpData : RawExpData
            Raw exposure data.
        """

        rawExpDirList = rawExpData.getRawExpDir()
        for rawExpDir in rawExpDirList:
            rawImgFiles = os.path.join(rawExpDir, "*.fits")
            self.wepCntlr.dataCollector.ingestImages(rawImgFiles)


if __name__ == "__main__":
    pass
