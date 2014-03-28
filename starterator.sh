#!/bin/bash
abspath="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
echo $abspath
starterator_path=`dirname $abspath`
echo $starterator_path

cd $starterator_path
<<<<<<< HEAD
git pull
<<<<<<< HEAD
=======
bzr pull
>>>>>>> 17a3d1a2438d87df8c80286abfafb97e14028e9f
$starterator_path/Starterator
=======
python $starterator_path/Starterator
>>>>>>> refactor
