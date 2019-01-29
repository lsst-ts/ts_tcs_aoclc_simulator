import os
import re
import numpy as np

from lsst.ts.phosim.Utility import getModulePath
from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.TeleFacade import TeleFacade


class WepPhosimCmpt(object):

    def __init__(self, phosimDir):
        """Initialization of WEP PhoSim component class.

        WEP: wavefront estimation pipeline.

        Parameters
        ----------
        phosimDir : str
            PhoSim directory.
        """

        self.metr = OpdMetrology()
        self.tele = TeleFacade()

        self.outputDir = ""
        self.outputImgDir = ""
        self.seedNum = None

        self.phosimParam = {"numPro": 1,
                            "e2ADC": 1}

        self._config(phosimDir)

    def _config(self, phosimDir):
        """Do the configuration of WEP PhoSim component class.

        Parameters
        ----------
        phosimDir : str
            PhoSim directory.
        """

        # Configurate the TeleFacade class
        configDataPath = self._getConfigDataPath()
        teleConfigFilePath = os.path.join(configDataPath, "telescopeConfig",
                                          "GT.inst")
        self.tele.setConfigFile(teleConfigFilePath)

        camDataDir = os.path.join(configDataPath, "camera")
        M1M3dataDir = os.path.join(configDataPath, "M1M3")
        M2dataDir = os.path.join(configDataPath, "M2")
        self.tele.setSubSysConfigDir(camDataDir=camDataDir,
                                     M1M3dataDir=M1M3dataDir,
                                     M2dataDir=M2dataDir,
                                     phosimDir=phosimDir)

        self.tele.setSensorOn(sciSensorOn=True, wfSensorOn=False,
                              guidSensorOn=False)
        self.tele.setInstName("lsst15")

    def _getConfigDataPath(self):
        """Get the configuration data path.

        Returns
        -------
        str
            Configuration data path
        """

        configDataPath = os.path.join(getModulePath(), "configData")

        return configDataPath

    def setOutputDir(self, outputDir):
        """Set the output directory.

        The output directory will be constructed if there is no existed one.

        Parameters
        ----------
        outputDir : str
            Output directory.
        """

        self._makeDir(outputDir)
        self.outputDir = outputDir

    def _makeDir(self, newDir, exist_ok=True):
        """Make the new directory.

        Super-mkdir; create a leaf directory and all intermediate ones. Works
        like mkdir, except that any intermediate path segment (not just the
        rightmost) will be created if it does not exist.

        Parameters
        ----------
        newDir : str
            New directory.
        exist_ok : bool, optional
            If the target directory already exists, raise an OSError if
            exist_ok is False. Otherwise no exception is raised. (the default
            is True.)
        """

        os.makedirs(newDir, exist_ok=exist_ok)

    def getOutputDir(self):
        """Get the output directory.

        Returns
        -------
        str
            Output directory.
        """

        return self.outputDir

    def setOutputImgDir(self, outputImgDir):
        """Set the output image directory.

        The output image directory will be constructed if there is no existed
        one.

        Parameters
        ----------
        outputImgDir : str
            Output image directory
        """

        self._makeDir(outputImgDir)
        self.outputImgDir = outputImgDir

    def getOutputImgDir(self):
        """Get the output image directory.

        Returns
        -------
        str
            Output image directory
        """

        return self.outputImgDir

    def setSeedNum(self, seedNum):
        """Set the seed number for the M1M3 mirror surface purturbation.

        Parameters
        ----------
        seedNum : int
            Seed number.
        """

        self.seedNum = int(seedNum)

    def getSeedNum(self):
        """Get the seed number for the M1M3 random surface purturbation.

        Returns
        -------
        int or None
            Seed number. None means there is no random purturbation.
        """

        return self.seedNum

    def setPhosimParam(self, numPro=1, e2ADC=1):
        """Set the PhoSim simulation parameters.

        Parameters
        ----------
        numPro : int, optional
            Number of processors. The value should be greater than 1. (the
            default is 1.)
        e2ADC : int, optional
            Whether to generate amplifier images (1 = true, 0 = false) (the
            default is 1.)
        """

        if (numPro > 0):
            self.phosimParam["numPro"] = int(numPro)

        if e2ADC in (0, 1):
            self.phosimParam["e2ADC"] = int(e2ADC)

    def getPhosimParam(self):
        """Get the PhoSim simulation parameters.

        Returns
        -------
        dict
            PhoSim simulation parameters.
        """

        return self.phosimParam

    def setSurveyParam(self, obsId=None, filterType=None, boresight=None,
                       zAngleInDeg=None, rotAngInDeg=None, mjd=None):
        """Set the survey parameters.

        Parameters
        ----------
        obsId : int, optional
            Observation Id. (the default is None.)
        filterType : FilterType, optional
            Active filter type. (the default is None.)
        boresight : tuple, optional
            Telescope boresight in (ra, decl). (the default is None.)
        zAngleInDeg : float, optional
            Zenith angle in degree. (the default is None.)
        rotAngInDeg : float, optional
            Camera rotation angle in degree between -90 and 90 degrees. (the
            default is None.)
        mjd : float, optional
            MJD of observation. (the default is None.)
        """

        self.tele.setSurveyParam(obsId=obsId, filterType=filterType,
                                 boresight=boresight, zAngleInDeg=zAngleInDeg,
                                 rotAngInDeg=rotAngInDeg, mjd=mjd)

    def addOpdFieldXYbyDeg(self, fieldXInDegree, fieldYInDegree):
        """Add the OPD new field X, Y in degree.

        OPD: optical path difference.

        Parameters
        ----------
        fieldXInDegree : float, list, or numpy.ndarray
            New field X in degree.
        fieldYInDegree : float, list, or numpy.ndarray
            New field Y in degree.
        """

        self.metr.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

    def accDofInUm(self, dofInUm):
        """Accumulate the aggregated degree of freedom (DOF) in um.

        idx 0-4: M2 dz, dx, dy, rx, ry
        idx 5-9: Cam dz, dx, dy, rx, ry
        idx 10-29: M1M3 20 bending modes
        idx 30-49: M2 20 bending modes

        Parameters
        ----------
        dofInUm : list or numpy.ndarray
            DOF in um.
        """

        self.tele.accDofInUm(dofInUm)

    def runPhoSim(self, argString):
        """Run the PhoSim program.

        Parameters
        ----------
        argString : str 
            Arguments for PhoSim.
        """

        self.tele.runPhoSim(argString)

    def getOpdArgsAndFilesForPhoSim(self, cmdFileName="opd.cmd",
                                    instFileName="opd.inst",
                                    logFileName="opdPhoSim.log",
                                    cmdSettingFileName="opdDefault.cmd",
                                    instSettingFileName="opdDefault.inst"):
        """Get the OPD calculation arguments and files for the PhoSim
        calculation.

        OPD: optical path difference.

        Parameters
        ----------
        cmdFileName : str, optional
            Physical command file name. (the default is "opd.cmd".)
        instFileName : str, optional
            OPD instance file name. (the default is "opd.inst".)
        logFileName : str, optional
            Log file name. (the default is "opdPhoSim.log".)
        cmdSettingFileName : str, optional
            Physical command setting file name. (the default is
            "opdDefault.cmd".)
        instSettingFileName : str, optional
            Instance setting file name. (the default is "opdDefault.inst".)

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        # Write the command file
        cmdFilePath = self._writePertAndCmdFiles(cmdSettingFileName,
                                                 cmdFileName)

        # Write the instance file
        instSettingFile = self._getInstSettingFilePath(instSettingFileName)
        instFilePath = self.tele.writeOpdInstFile(
            self.outputDir, self.metr, instSettingFile=instSettingFile,
            instFileName=instFileName)

        # Get the argument to run the PhoSim
        argString = self._getPhoSimArgs(logFileName, instFilePath, cmdFilePath)

        return argString

    def _writePertAndCmdFiles(self, cmdSettingFileName, cmdFileName):
        """Write the physical perturbation and command files.

        Parameters
        ----------
        cmdSettingFileName : str
            Physical command setting file name.
        cmdFileName : str
            Physical command file name.

        Returns
        -------
        str
            Command file path.
        """

        # Write the perturbation file
        pertCmdFilePath = self.tele.writePertBaseOnConfigFile(
            self.outputDir, seedNum=self.seedNum, saveResMapFig=True,
            pertCmdFileName="pert.cmd")

        # Write the physical command file
        configDataPath = self._getConfigDataPath()
        cmdSettingFile = os.path.join(configDataPath, "cmdFile",
                                      cmdSettingFileName)
        cmdFilePath = self.tele.writeCmdFile(
            self.outputDir, cmdSettingFile=cmdSettingFile,
            pertFilePath=pertCmdFilePath, cmdFileName=cmdFileName)

        return cmdFilePath

    def _getInstSettingFilePath(self, instSettingFileName):
        """Get the instance setting file path.

        Parameters
        ----------
        instSettingFileName : str
            Instance setting file name.

        Returns
        -------
        str
            Instance setting file path.
        """

        configDataPath = self._getConfigDataPath()
        instSettingFile = os.path.join(configDataPath, "instFile",
                                       instSettingFileName)

        return instSettingFile

    def _getPhoSimArgs(self, logFileName, instFilePath, cmdFilePath):
        """Get the arguments needed to run the PhoSim.

        Parameters
        ----------
        logFileName : str
            Log file name.
        instFilePath: str
            Instance file path.
        cmdFilePath : str
            Physical command file path.

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        numPro = self.phosimParam["numPro"]
        e2ADC = self.phosimParam["e2ADC"]
        logFilePath = os.path.join(self.outputImgDir, logFileName)

        argString = self.tele.getPhoSimArgs(
            instFilePath, extraCommandFile=cmdFilePath, numPro=numPro,
            outputDir=self.outputImgDir, e2ADC=e2ADC, logFilePath=logFilePath)

        return argString

    def getComCamOpdArgsAndFilesForPhoSim(
            self, cmdFileName="opd.cmd", instFileName="opd.inst",
            logFileName="opdPhoSim.log", cmdSettingFileName="opdDefault.cmd",
            instSettingFileName="opdDefault.inst"):
        """Get the OPD calculation arguments and files of ComCam for the PhoSim
        calculation.

        OPD: optical path difference.
        ComCam: commissioning camera.

        Parameters
        ----------
        cmdFileName : str, optional
            Physical command file name. (the default is "opd.cmd".)
        instFileName : str, optional
            OPD instance file name. (the default is "opd.inst".)
        logFileName : str, optional
            Log file name. (the default is "opdPhoSim.log".)
        cmdSettingFileName : str, optional
            Physical command setting file name. (the default is
            "opdDefault.cmd".)
        instSettingFileName : str, optional
            Instance setting file name. (the default is "opdDefault.inst".)

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        # Set the default ComCam OPD field positions
        self.metr.setDefaultComcamGQ()

        argString = self.getOpdArgsAndFilesForPhoSim(
            cmdFileName=cmdFileName, instFileName=instFileName,
            logFileName=logFileName, cmdSettingFileName=cmdSettingFileName,
            instSettingFileName=instSettingFileName)

        return argString

    def getStarArgsAndFilesForPhoSim(self, skySim,
                                     cmdFileName="star.cmd",
                                     instFileName="star.inst",
                                     logFileName="starPhoSim.log",
                                     simSeed=1000,
                                     cmdSettingFileName="starDefault.cmd",
                                     instSettingFileName="starSingleExp.inst"):
        """Get the star calculation arguments and files for the PhoSim
        calculation.

        Parameters
        ----------
        skySim : SkySim
            Sky simulator
        cmdFileName : str, optional
            Physical command file name. (the default is "star.cmd".)
        instFileName : str, optional
            Star instance file name. (the default is "star.inst".)
        logFileName : str, optional
            Log file name. (the default is "starPhoSim.log".)
        simSeed : int, optional
            Random number seed. (the default is 1000)
        cmdSettingFileName : str, optional
            Physical command setting file name. (the default is
            "starDefault.cmd".)
        instSettingFileName : str, optional
            Instance setting file name. (the default is "starSingleExp.inst".)

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        # Write the command file
        cmdFilePath = self._writePertAndCmdFiles(cmdSettingFileName,
                                                 cmdFileName)

        # Write the instance file
        instSettingFile = self._getInstSettingFilePath(instSettingFileName)
        instFilePath = self.tele.writeStarInstFile(self.outputDir, skySim,
                          simSeed=simSeed, sedName="sed_flat.txt",
                          instSettingFile=instSettingFile,
                          instFileName=instFileName)

        # Get the argument to run the PhoSim
        argString = self._getPhoSimArgs(logFileName, instFilePath, cmdFilePath)

        return argString

    def calcOpdPssn(self):
        """Calculate the PSSN of OPD.

        OPD: optical path difference.
        PSSN: normalized point source sensitivity.

        Returns
        -------
        list
            PSSN list.
        """

        opdFileList = self._getOpdFileInDir(self.outputImgDir)

        wavelengthInUm = self.tele.WAVELENGTH_IN_NM * 1e-3
        pssnList = []
        for opdFile in opdFileList:
            pssn = self.metr.calcPSSN(wavelengthInUm, opdFitsFile=opdFile)
            pssnList.append(pssn)

        return pssnList

    def _getOpdFileInDir(self, opdDir):
        """Get the OPD files in the directory.

        OPD: optical path difference.

        Parameters
        ----------
        opdDir : str
            OPD file directory.

        Returns
        -------
        list
            List of OPD files.
        """

        opdFileList = []
        fileList = self._getFileInDir(opdDir)
        for file in fileList:
            fileName = os.path.basename(file)
            m = re.match(r"\Aopd_\d+_(\d+).fits.gz", fileName)
            if (m is not None):
                opdFileList.append(file)

        return opdFileList

    def _getFileInDir(self, fileDir):
        """Get the files in the directory.

        Parameters
        ----------
        fileDir : str
            File directory.

        Returns
        -------
        list
            List of files.
        """

        fileList = []
        for name in os.listdir(fileDir):
            filePath = os.path.join(fileDir, name)
            if os.path.isfile(filePath):
                fileList.append(filePath)

        return fileList

    def calcComCamGQeffFwhm(self, pssnList):
        """Calculate the GQ effective FWHM of ComCam.

        GQ: Gaussian quadrature.
        FWHM: full width and half maximum.
        PSSN: normalized point source sensitivity.
        ComCam: commissioning camera.

        Parameters
        ----------
        pssnList : list
            List of PSSN.

        Returns
        -------
        float
            GQ effective FWHM of ComCam.
        """

        # Calculate the list of effective FWHM
        effFwhmList = []
        for pssn in pssnList:
            effFwhm = self.metr.calcFWHMeff(pssn)
            effFwhmList.append(effFwhm)

        # Calculate the GQ effectice FWHM
        comcamWtRatio = np.ones(9)
        self.metr.setWeightingRatio(comcamWtRatio)
        gqEffFwhm = self.metr.calcGQvalue(effFwhmList)

        return gqEffFwhm


if __name__ == "__main__":
    pass
