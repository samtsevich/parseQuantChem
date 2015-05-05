#!/bin/bash
for fullpath in *
do
    filename="${fullpath##*/}"                      # Strip longest match of */ from start
    dir="${fullpath:0:${#fullpath} - ${#filename}}" # Substring from 0 thru pos of filename
    base="${filename%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
    ext="${filename:${#base} + 1}"                  # Substring from len of base thru end

    if [ -d "$fullpath" ]; then

        # echo $fullpath >> 1.txt

        for path in $fullpath/*.out
        do
            filename1="${path##*/}"                      # Strip longest match of */ from start
            base1="${path%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
            echo $path >> 1.txt
            cp -r $path ~/res/$filename1
        done

        # for path in $fullpath/*.gjf
        # do
        #     filename1="${path##*/}"                      # Strip longest match of */ from start
        #     base1="${path%.[^.]*}"                       # Strip shortest match of . plus at least one non-dot char from end
        #     echo $path >> 1.txt
        #     mv $path ~/res/$filename1
        # done

        echo "_____________" >> 1.txt
    fi
done

echo ended script
exit 0
