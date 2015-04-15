#!/bin/bash
declare -i i
i=1
echo "    " > "1.txt"
for fullpath in *.xyz
do
    filename="${fullpath##*/}"                      # Strip longest match of */ from start
    dir="${fullpath:0:${#fullpath} - ${#filename}}" # Substring from 0 thru pos of filename
    base="${filename%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
    ext="${filename:${#base} + 1}"                  # Substring from len of base thru end

    mv $filename "Ca2HPO4_$i.xyz"
    echo "$filename     Ca2HPO4_$i" >> "1.txt"
    i=i+1
done

echo ended script
exit 0
