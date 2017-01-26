#!/bin/bash
##
## wrapper script for bwa - just wrap in a time command
##

USE_PATH_TO_PICARD=${PATH_TO_PICARD:-"/usr/picard/"}
time java -Xmx4g -jar "${USE_PATH_TO_PICARD}/picard.jar" "$@"
