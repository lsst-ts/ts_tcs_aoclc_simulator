#!/usr/bin/env groovy

pipeline {

    agent {
        // Use the docker to assign the Python version.
        // Use the label to assign the node to run the test.
        // The nodes in T&S teams is 'jenkins-el7-1'.
        // It is recommended by SQUARE team do not add the label.
        docker {
            image 'lsstts/aos_aoclc:w_2019_02'
            args '-u root'
        }
    }

    triggers {
        pollSCM('H * * * *')
    }

    environment {
        //Declare the aos module path
        AOS_MODULE_PATH="/home/lsst/aos_repos"
        // Use the double quote instead of single quote
        // Add the PYTHONPATH
        PYTHONPATH="${env.WORKSPACE}/python:${AOS_MODULE_PATH}/ts_tcs_wep/python:${AOS_MODULE_PATH}/ts_tcs_ofcPython/python:${AOS_MODULE_PATH}/ts_tcs_wep_phosim/python"
        // XML report path
        XML_REPORT="jenkinsReport/report.xml"
        // Module name used in the pytest coverage analysis
        MODULE_NAME="lsst.ts.aoclcSim"
    }

    stages {

        stage('Unit Tests and Coverage Analysis') { 
            steps {
                // Direct the HOME to WORKSPACE for pip to get the
                // installed library.
                // 'PATH' can only be updated in a single shell block.
                // We can not update PATH in 'environment' block.
                // Pytest needs to export the junit report. 
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        source /opt/rh/devtoolset-6/enable
                        source /opt/lsst/loadLSST.bash
                        setup sims_catUtils -t sims_w_2019_02
                        cd ${AOS_MODULE_PATH}/obs_lsst
                        setup -k -r .
                        cd ${AOS_MODULE_PATH}/phosim_utils
                        setup -k -r .
                        cd ${env.WORKSPACE}
                        pytest --cov-report html --cov=${env.MODULE_NAME} --junitxml=${env.WORKSPACE}/${env.XML_REPORT} ${env.WORKSPACE}/tests/*.py
                    """
                }
            }
        }
    }

    post {        
        always {
            // The path of xml needed by JUnit is relative to
            // the workspace.
            junit 'jenkinsReport/*.xml'

            // Publish the HTML report
            publishHTML (target: [
                allowMissing: false,
                alwaysLinkToLastBuild: false,
                keepAll: true,
                reportDir: 'htmlcov',
                reportFiles: 'index.html',
                reportName: "Coverage Report"
              ])
        }

        cleanup {
            // clean up the workspace
            deleteDir()
        }  
    }
}
