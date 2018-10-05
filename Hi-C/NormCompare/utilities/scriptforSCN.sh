#!/bin/bash

usage() { echo "Usage: $0 [-f <string>]" 1>&2; exit 1; }

while getopts ":f:" opt; do
    case "${opt}" in
        f)
            f=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${f}" ]; then
    usage
fi

matlab -nodesktop -nosplash -nodisplay -r "try, run('${f}'); catch err, disp(err.message); exit(1); end; exit(0)"
