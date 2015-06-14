#!/bin/bash

octave -q runOctave.m

cp -r resCoef ../../resCoef 
cp -r resEner ../../resEner 

echo ended script
exit 0
