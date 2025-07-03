#!/bin/bash
#set -x
idfile="${1}"
while read line ; do
    #echo "${line}" | sed -E 's/^.type":"item","id":\("Q\).*$/\1/'
    theid=`echo "${line}" | sed -n 's/^."type":"item","id":\("Q[0-9]*"\).*$/\1/p'`
    if ! test -z $theid ; then 
        if grep "${theid}" "${idfile}" >/dev/null ; then
            echo "${line}"
        fi
    fi
done
