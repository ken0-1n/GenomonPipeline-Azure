#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Fastq_splitter(Stage_task):

    task_name = "fastq_splitter"

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


to_val=`ls {target_dir}/*_{pair_num}{ext} | wc -l`

echo $to_val

input_files=""
for i in `seq 1 ${{to_val}}`; do
    input_files="${{input_files}} {target_dir}/${{i}}_{pair_num}{ext}"
done    

if [ "_{fastq_filter}" = "_True" ]; then

    if [ "_{ext}" = "_.gz" ]; then
        gzip -dc $input_files | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' > {target_dir}/{pair_num}.in.fastq || exit $?
    else
        cat $input_files | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' > {target_dir}/{pair_num}.in.fastq || exit $?
    fi
else
    if [ "_{ext}" = "_.gz" ]; then
        gzip -dc $input_files > {target_dir}/{pair_num}.in.fastq || exit $?
    else
        cp $input_files {target_dir}/{pair_num}.in.fastq || exit $?
    fi
fi


"""

    def __init__(self, qsub_option, script_dir):
        super(Fastq_splitter, self).__init__(qsub_option, script_dir)


