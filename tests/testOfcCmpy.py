import os
import numpy as np
import unittest

from lsst.ts.wep.Utility import FilterType
from lsst.ts.ofc.IterDataReader import IterDataReader

from lsst.ts.aoclcSim.Utility import getModulePath
from lsst.ts.aoclcSim.OfcCmpt import OfcCmpt


class TestOfcCmpt(unittest.TestCase):
    """ Test the OfcCmpt class."""

    def setUp(self):

        # Set the OfcCmpt object
        self.ofcCmpt = OfcCmpt()
        self.ofcCmpt.setFilter(FilterType.REF)
        self.ofcCmpt.setRotAng(0)

        # Set the interation data reader
        iterDataDir = os.path.join(getModulePath(), "tests", "testData",
                                   "comcamIterData")
        self.iterDataReader = IterDataReader(iterDataDir)

    def testCalcAggDofForPhoSim(self):

        sensorNameList = self._getComCamSensorNameList()
        pssnSensorNameList = sensorNameList

        maxIterNum = 5
        for iterNum in range(0, maxIterNum):

            # Get the PSSN. The final one value in data is the GQ PSSN.
            # Therefore, only take the first 9 values.
            pssn = self.iterDataReader.getPssn(iterNum)
            pssn = pssn[0:9]

            # Set the gain value based on PSSN
            self.ofcCmpt.setGainByPSSN(pssn, sensorNameList)

            # Get the wavefront error
            wfErr = self.iterDataReader.getWfsErr(iterNum)

            # Calculate the DOF
            dof = self.ofcCmpt.calcAggDofForPhoSim(wfErr, sensorNameList)

            # Get the answer of DOF from the test data
            dofAns = self.iterDataReader.getDof(iterNum + 1)

            # Calculate the difference and do the assertion
            delta = np.sum(np.abs(dof - dofAns))
            self.assertLess(delta, 0.002)

    def _getComCamSensorNameList(self):

        sensorNameList = ["R22_S00", "R22_S01", "R22_S02", "R22_S10", "R22_S11",
                          "R22_S12", "R22_S20", "R22_S21", "R22_S22"]
        return sensorNameList


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
