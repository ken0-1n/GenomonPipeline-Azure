#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Fusionfusion(Stage_task):

    task_name = "fusionfusion"

    script_template = """
#!/bin/bash
#
# Set SGE
#
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

# set python environment
if [ "_{pythonhome}" != "_" ];then
    export PYTHONHOME={pythonhome}
    export PATH=$PYTHONHOME/bin:$PATH
fi
if [ "_{ld_library_path}" != "_" ];then
    export LD_LIBRARY_PATH={ld_library_path}
fi
if [ "_{pythonpath}" != "_" ];then
    export PYTHONPATH={pythonpath}
fi

{fusionfusion} --star {chimeric_sam} --out {output_prefix} --param {param_file}
"""

    def __init__(self, qsub_option, script_dir):
        super(Fusionfusion, self).__init__(qsub_option, script_dir)
