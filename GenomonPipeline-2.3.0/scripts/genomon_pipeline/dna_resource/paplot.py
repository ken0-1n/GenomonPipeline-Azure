#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_PA_Plot(Stage_task):

    task_name = "paplot"

    script_template = """
#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv
# set python environment
if [ "_{pythonhome}" != "_" ];then
    export PYTHONHOME={pythonhome}
    export PATH=$PYTHONHOME/bin:$PATH
fi
if [ "_{pythonpath}" != "_" ];then
    export PYTHONPATH={pythonpath}
fi

if test "{inputs_qc}" != ""; then
  {pa_plot} qc "{inputs_qc}" {output_dir} {title} --config_file {config_file} --remarks "{remarks}"
fi
if test "{inputs_sv}" != ""; then
  {pa_plot} sv "{inputs_sv}" {output_dir} {title} --config_file {config_file} --remarks "{remarks}"
fi
if test "{inputs_mutation}" != ""; then
  {pa_plot} mutation "{inputs_mutation}" {output_dir} {title} --config_file {config_file} --remarks "{remarks}"
fi
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PA_Plot, self).__init__(qsub_option, script_dir)

        
