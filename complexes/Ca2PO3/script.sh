#!/bin/bash
for fullpath in *.xyz
do
    filename="${fullpath##*/}"                      # Strip longest match of */ from start
    dir="${fullpath:0:${#fullpath} - ${#filename}}" # Substring from 0 thru pos of filename
    base="${filename%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
    ext="${filename:${#base} + 1}"                  # Substring from len of base thru end

    echo $base
    echo $filename

	# sed -i '1s/^/ 7/' $fullpath
done

echo ended script
exit 0
