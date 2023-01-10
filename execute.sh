#!/bin/bash
# shell script to download dataset from terra and execute the pipeline

# docker build -t amir_docker -f Docker .

# ###############
# Download data
# ###############
set +eu && PS1=dummy && . /opt/conda/etc/profile.d/conda.sh; conda activate Seurat_conda; gsutil -m cp -r gs://terra-featured-workspaces/Cumulus/cellranger_output/ /temp

# ##########
# process input data into json format
# ###########

for i in $(find "$(cd /temp; pwd)"  -name "*.h5"); do dir=$(dirname $i); cp $i /data/${dir##*/}.h5; done

echo -n '{ "pipeline.h5_files": ' > /data/input.json
find "$(cd /data; pwd)" -name "*.h5" | jq -R -s -c 'split("\n")[:-1]' >> /data/input.json
echo "}" >> /data/input.json


# #############
# Run the workflow
# ##############

cd /temp; cromwell run -i /data/input.json --options  /script/Genentech_Assignment/cromwell.conf /script/Genentech_Assignment/pipeline.wdl





