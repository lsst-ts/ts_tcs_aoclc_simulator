# Active Optics Closed-Loop Control (AOCLC) Simulator

*This module simulates the active optics closed-loop control (aoclc) with PhoSim.*

## 1. Version History

The version history is [here](./doc/VersionHistory.md).

*Author: Te-Wei Tsai* \
*Date: 2-28-2019*

## 2. Platform

- *CentOS 7*
- *python: 3.6.6*
- *scientific pipeline (newinstall.sh from master branch)*
- *phosim_syseng4 (branch: aos, tag: firstdonuts)*

## 3. Needed Package

- *lsst_sims (tag: sims_w_2019_02)*
- *lsst_distrib (tag: w_2019_02)*
- *obs_lsst - master branch (commit: 69b4a98)*
- *phosim_utils - master branch (commit: b8d87d9)*
- *ts_tcs_wep - master branch (commit: 00a021b)*
- *ts_tcs_wep_phosim - master branch (commit: 6e4d997)*
- *ts_tcs_ofcPython - master branch (commit: 53ee625)*

## 4. Install the LSST Packages, obs_lsst, and phosim_utils

*1. Setup the LSST environment by `source $LSST_DIR/loadLSST.bash`. LSST_DIR is the directory of scientific pipeline.* \
*2. Install the lsst_sims by `eups distrib install lsst_sims -t sims_w_2019_02`.* \
*3. Install the lsst_distrib by `eups distrib install lsst_distrib -t w_2019_02`.* \
*4. Fix the path by `curl -sSL https://raw.githubusercontent.com/lsst/shebangtron/master/shebangtron | python`. The [shebangtron repo](https://github.com/lsst/shebangtron) has the further discussion of this.* \
*5. Clone the repository of [obs_lsst](https://github.com/lsst/obs_lsst) to some other directory. Under the obs_lsst directory, use `setup -k -r .` to setup the package in eups and use `scons` to build the module. It is noted that the build process is only needed for the first time.* \
*6. Do the step 5 for the repository of [phosim_utils](https://github.com/lsst-dm/phosim_utils.git).*

## 5. Pull the Built Image from Docker Hub

*Pull the built docker image by `docker pull lsstts/aos_aoclc:w_2019_02`. The scientific pipeline and aos-related packages are installed already. For the details of docker image, please follow the [docker aoclc image](https://hub.docker.com/r/lsstts/aos_aoclc).*

## 6. Use of Module

*1. Setup the DM environment:*

```bash
source $path_of_lsst_scientific_pipeline/loadLSST.bash
setup sims_catUtils -t $user_defined_tag -t sims_w_2019_02
```

*2. Setup the WEP environment:* \
`export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_wep/python`

*3. Setup the OFC environment:* \
`export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_ofcPython/python`

*4. Setup the wepPhoSim environment:* \
`export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_wep_phosim/python`

*5. Setup the AOCLC simulator environment:* \
`export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_aoclc_simulator/python`

*6. Setup the PhoSim path variable:* \
`export PHOSIMPATH=$path_to_phosim_directory`

*7. If use the docker image, only need to do the steps 5 and 6.*

## 7. Content

*This module contains the following classes and functions ([class diagram](./doc/aoclcClassDiag.png)):*

- **WepPhosimCmpt**: High-level component to use the module of ts_tcs_wep_phosim.
- **WepCmpt**: High-level component to use the module of ts_tcs_wep.
- **OfcCmpt**: High-level component to use the module of ts_tcs_ofcPython.

## 8. Example Script

- **opdClosedLoop.py**: Closed-loop simulation in the optical path difference (OPD) level, which means only the classes of WepPhosimCmpt and OfcCmpt are used.
- **comcamClosedLoop.py**: Closed-loop simulation of commionning camera. There are 9 stars on the center of each CCD.