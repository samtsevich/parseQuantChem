#!/bin/bash

for fullpath in *.out
do
    filename="${fullpath##*/}"                      # Strip longest match of */ from start
    dir="${fullpath:0:${#fullpath} - ${#filename}}" # Substring from 0 thru pos of filename
    base="${filename%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
    ext="${filename:${#base} + 1}"                  # Substring from len of base thru end

    mv $filename "../resCharges/$filename"
done

echo ended script
exit 0