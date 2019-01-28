# Active Optics Closed-Loop Control (AOCLC) Simulator

*This module simulates the active optics closed-loop control (aoclc) with PhoSim.*

## 1. Version History

*Version 0.1*
<br/>
*Simulated the aoclc with the commionning camera (ComCam).*
<br/>

*Author: Te-Wei Tsai*
<br/>
*Date: 1-28-2019*

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
- *ts_tcs_wep - develop branch (commit: d59002a)*
- *ts_tcs_wep_phosim - develop branch (commit: 6e4d997)*
- *ts_tcs_ofcPython - develop branch (commit: 107811e)*

## 4. Pull the Built Image from Docker Hub

*Pull the built docker image by `docker pull lsstts/aos_aoclc:w_2019_02`. The scientific pipeline and aos-related packages are installed already. For the details of docker image, please follow the [docker aoclc image](https://hub.docker.com/r/lsstts/aos_aoclc).*