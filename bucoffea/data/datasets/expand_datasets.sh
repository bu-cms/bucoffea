#!/bin/bash


expand(){
    INFILE=${1}
    COND=${2}
    for STUB in $(cat ${INFILE}); do
        dasgoclient --query="dataset=/${STUB}/${COND}/NANOAOD*"
    done
}

INFILE="datasetnames.txt"
expand "${INFILE}" "RunIIFall17*1June*" | tee datasets_2017.txt
expand "${INFILE}" "RunIIAutumn18*1June*" | tee datasets_2018.txt