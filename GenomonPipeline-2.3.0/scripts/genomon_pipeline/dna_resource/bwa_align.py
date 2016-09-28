#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Bwa_align(Stage_task):

    task_name = "bwa_align"

    script_template = """
#!/usr/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

{bwa} mem {bwa_params} {ref_fa} {input_dir}/1.in.fastq {input_dir}/2.in.fastq > {output_dir}/{sample_name}.bwa.sam
status=("${{PIPESTATUS[@]}}")
if [ ${{status[0]}} -ne 0 ]; then
    exit ${{status[0]}}
fi
if [ ${{status[1]}} -ne 0 ]; then
    exit ${{status[1]}}
fi

{biobambam}/bamsort index=1 level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindentonly=1 calmdnmreference={ref_fa} tmpfile={output_dir}/{sample_name}.sorted.bam.tmp inputformat=sam indexfilename={output_dir}/{sample_name}.sorted.bam.bai I={output_dir}/{sample_name}.bwa.sam O={output_dir}/{sample_name}.sorted.bam

"""

    def __init__(self, qsub_option, script_dir):
        super(Bwa_align, self).__init__(qsub_option, script_dir)


