import os
import lsst.ts.aoclcSim


def getModulePath(module=lsst.ts.aoclcSim, startIdx=1, endIdx=-4):
    """Get the path of module.

    Parameters
    ----------
    module : str, optional
        Module name. (the default is lsst.ts.aoclcSim.)
    startIdx : int, optional
        Start index. (the default is 1.)
    endIdx : int, optional
        End index. (the default is -4.)

    Returns
    -------
    str
        Directory path of module based on the start and end indexes.
    """

    # Get the path of module
    modulePathList = os.path.dirname(module.__file__).split(
                                os.sep)[int(startIdx):int(endIdx)]
    modulePath = os.path.join(os.sep, *modulePathList)

    return modulePath


def getPhoSimPath(phosimPathVar="PHOSIMPATH"):
    """Ge the PhoSim path from the environment variables.

    Parameters
    ----------
    phosimPathVar : str, optional
        PhoSim path variable name. (the default is "PHOSIMPATH".)

    Returns
    -------
    str
        PhoSim path.

    Raises
    ------
    ValueError
        Please set the 'PHOSIMPATH' environment variable.
    """

    phosimPath = None
    try:
        phosimPath = os.environ[phosimPathVar]
    except Exception as KeyError:
        raise ValueError("Please set the '%s' environment variable."
                         % phosimPathVar)

    return phosimPath


if __name__ == "__main__":
    pass
