#!/bin/bash

#read flags, use defaults

export PATH=/opt/orca:$PATH
export RSH_COMMAND="/usr/bir/ssh -x"
INPUT=`pwd`

#set long flags
for arg in "$@"; do
    shift
    case "$arg" in
	"--input") set -- "$@" "-i" ;;
	*)        set -- "$@" "$arg"
        esac
done

#set internal flags
while getopts ":i:" Option
do
  case $Option in
    i) INPUT=`pwd`/${OPTARG}
       cd $INPUT ;;
    esac
    done
    shift $(($OPTIND - 1))

exec 1>$INPUT/log.log
exec 2>$INPUT/err_log.log

#start computation
for filename in `find $INPUT -type f -name \*.inp`
do
	/home/rudenko/orca/orca $filename > $(basename "$filename" | cut -d. -f1).out
done