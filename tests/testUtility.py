import os
import unittest

from lsst.ts.aoclcSim.Utility import getPhoSimPath


class TestUtility(unittest.TestCase):
    """ Test the functions in Utility."""

    def testGetPhoSimPathNotExist(self):

        self.assertRaises(ValueError, getPhoSimPath, "WRONGPATH")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
