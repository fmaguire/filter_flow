#!/bin/bash


pwd=`pwd`
for i in `ls conda_envs/*.yml`; do
   fileName=`basename $i`;
   envName=NXF_FILTERING_"${fileName%.*}"
   envPath=$pwd/conda_envs/$envName
   #echo $i;
   #echo $envName;
   if ! { conda env list | grep $envName; } >/dev/null 2>&1; then
	  conda create -p $envPath -y;
	  conda env update -p $envPath --file $pwd/$i;
   else
	  echo "Env $envName already exists, passing"
   fi
done;