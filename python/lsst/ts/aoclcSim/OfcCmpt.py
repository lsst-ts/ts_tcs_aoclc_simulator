import os
import numpy as np

from lsst.ts.wep.Utility import FilterType

from lsst.ts.ofc.Utility import getModulePath, InstName, DofGroup
from lsst.ts.ofc.DataShare import DataShare
from lsst.ts.ofc.OptStateEstiDataDecorator import OptStateEstiDataDecorator
from lsst.ts.ofc.OptCtrlDataDecorator import OptCtrlDataDecorator
from lsst.ts.ofc.OptStateEsti import OptStateEsti
from lsst.ts.ofc.OptCtrl import OptCtrl
from lsst.ts.ofc.ZTAAC import ZTAAC
from lsst.ts.ofc.CamRot import CamRot


class OfcCmpt(object):

    DEFAULT_GAIN = 0.7
    FWHM_THRESHOLD_IN_ARCSEC = 0.2

    def __init__(self):
        """Initialization of OFC component class. At this moment, it only
        supports the ComCam as a minimum requirement."""

        self.ztaac = self._configZTAAC(InstName.COMCAM)
        self.camRot = CamRot(rotAngInDeg=0)

    def _configZTAAC(self, instName):
        """Configurate the ZTAAC.

        ZTAAC: Zernike to actuator adjustment calculator.

        Parameters
        ----------
        instName : InstName
            Instrument name.

        Returns
        -------
        ZTAAC
            configurated ZTAAC object.
        """

        # Prepare the data object for the ZTAAC to use
        dataShare = DataShare()
        configDataPath = self._getConfigDataPath()
        dataShare.config(configDataPath, instName=instName)

        optStateEstiData = OptStateEstiDataDecorator(dataShare)
        optStateEstiData.configOptStateEstiData()

        mixedData = OptCtrlDataDecorator(optStateEstiData)
        mixedData.configOptCtrlData(configFileName="optiPSSN_x00.ctrl")

        # Instantiate the ZTAAC object with the configured objects.
        ztaac = ZTAAC(OptStateEsti(), OptCtrl(), mixedData)

        # Set the state 0 and state
        ztaac.setState0FromFile(state0InDofFileName="state0inDof.txt")
        ztaac.setStateToState0()

        # Configure the parameters
        ztaac.config(filterType=FilterType.REF, defaultGain=self.DEFAULT_GAIN,
                     fwhmThresholdInArcsec=self.FWHM_THRESHOLD_IN_ARCSEC)

        return ztaac

    def _getConfigDataPath(self):
        """Get the configuration data path.

        Returns
        -------
        str
            Configuration data path
        """

        configDataPath = os.path.join(getModulePath(), "configData")

        return configDataPath

    def setFilter(self, filterType):
        """Set the current filter.

        Parameters
        ----------
        filterType : FilterType
            The new filter configuration to use for OFC data processing.
        """

        self.ztaac.setFilter(filterType)

    def getFilter(self):
        """Get the current filter.

        Returns
        -------
        FilterType
            The current filter configuration to use for OFC data processing.
        """

        return self.ztaac.getFilter()

    def setRotAng(self, rotAngInDeg):
        """Set the camera rotation angle in degree.

        Parameters
        ----------
        rotAngInDeg : float
            The camera rotation angle in degree (-90 to 90).
        """

        self.camRot.setRotAng(rotAngInDeg)

    def getRotAng(self):
        """Get the camera rotation angle in degree.

        Returns
        -------
        float
            The camera rotation angle in degree.
        """

        return self.camRot.getRotAng()

    def setGainByPSSN(self, pssn, sensorNameList):
        """Set the gain value based on PSSN.

        PSSN: Normalized point spread function.

        Parameters
        ----------
        pssn : numpy.ndarray or list
            PSSN.
        sensorNameList : list[str]
            List of abbreviated sensor names.
        """

        self.ztaac.setGainByPSSN(pssn, sensorNameList)

    def setGain(self, gain):
        """Set the gain value.

        Parameters
        ----------
        gain : float
            Gain value in the feedback. It should be in the range of 0 and 1.
        """

        self.ztaac.setGain(gain)

    def getState0(self):
        """Get the state 0 in degree of freedom (DOF).

        Returns
        -------
        numpy.ndarray
            State 0 in DOF.
        """

        return self.ztaac.getState0()

    def calcAggDofForPhoSim(self, wfErr, sensorNameList):
        """Calculate the aggregated DOF for PhoSim simulation.

        DOF: Degree of freedom.

        Parameters
        ----------
        wfErr : numpy.ndarray
            Wavefront error.
        sensorNameList : list[str]
            List of abbreviated sensor names.

        Returns
        -------
        numpy.ndarray
            Aggregated DOF.
        """

        # Calculate the uk based on the control algorithm.
        uk = self.ztaac.estiUkWithGain(wfErr, sensorNameList)

        # Consider the camera rotation.
        rotUk = self.ztaac.rotUk(self.camRot, uk)

        # Aggregate the rotated uk.
        self.ztaac.aggState(rotUk)

        # Collect the DOF for the comparison
        dof = []
        for dofGroup in DofGroup:
            # Get the DOF for each group
            dofOfGroup = self.ztaac.getGroupDof(dofGroup)
            dof = np.append(dof, dofOfGroup)

        # In the simulation mode, the PhoSim needs the total state
        # of telescope to do the perturbation.
        dof += self.ztaac.getState0()

        return dof


if __name__ == "__main__":
    pass
