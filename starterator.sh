#!/bin/bash
abspath="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
echo $abspath
starterator_path=`dirname $abspath`
echo $starterator_path

cd $starterator_path

git pull
$starterator_path/Starterator
python $starterator_path/starterator/Starterator
