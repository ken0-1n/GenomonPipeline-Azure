#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class SV_filt(Stage_task):

    task_name = "sv_filt"

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
    export PATH={htslib}:$PYTHONHOME/bin:$PATH
fi
if [ "_{ld_library_path}" != "_" ];then
    export LD_LIBRARY_PATH={ld_library_path}
fi
if [ "_{pythonpath}" != "_" ];then
    export PYTHONPATH={pythonpath}
fi

{genomon_sv} filt {input_bam} {output_prefix} {reference_genome} {annotation_dir} {param} || exit $?

mv {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.txt.tmp || exit $?

echo -e "{meta_info}" > {output_prefix}.genomonSV.result.txt || exit $?

cat {output_prefix}.genomonSV.result.txt.tmp >> {output_prefix}.genomonSV.result.txt || exit $?

rm -rf {output_prefix}.genomonSV.result.txt.tmp
 
{sv_utils} filter {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.filt.txt.tmp {sv_utils_annotation_dir} {sv_utils_param} || exit $?

mv {output_prefix}.genomonSV.result.filt.txt.tmp {output_prefix}.genomonSV.result.filt.txt
"""

    def __init__(self, qsub_option, script_dir):
        super(SV_filt, self).__init__(qsub_option, script_dir)


