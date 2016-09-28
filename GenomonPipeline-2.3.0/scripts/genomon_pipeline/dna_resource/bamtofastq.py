#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Bam2Fastq(Stage_task):

    task_name = "bam2fastq"

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

bams=( `echo "{input_bam}" | tr -s ';' ' '`)

if [ ${{#bams[@]}} -eq 1 ]; then
    bam=${{bams[0]}}
    echo $bam
    {biobambam}/bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename=${{bam}} F={f1_name} F2={f2_name} T={t} S={s} O={o1_name} O2={o2_name} || exit $?    

else
    echo -n > {f1_name}
    echo -n > {f2_name}
    echo -n > {t}
    echo -n > {s}
    echo -n > {o1_name}
    echo -n > {o2_name}

    for bam in ${{bams[@]}}; do
        echo $bam
        {biobambam}/bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename=${{bam}} F={f1_name}.tmp F2={f2_name}.tmp T={t}.tmp S={s}.tmp O={o1_name}.tmp O2={o2_name}.tmp || exit $?
        cat {f1_name}.tmp >> {f1_name} || exit $?
        cat {f2_name}.tmp >> {f2_name} || exit $?
        if [ -s {t}.tmp ]; then
            cat {t}.tmp >> {t} || exit $?
        fi
        if [ -s {s}.tmp ]; then
            cat {s}.tmp >> {s} || exit $?
        fi
        if [ -s {o1_name}.tmp ]; then
            cat {o1_name}.tmp >> {o1_name} || exit $?
        fi
        if [ -s {o2_name}.tmp ]; then
            cat {o2_name}.tmp >> {o2_name} || exit $?
        fi
        rm {f1_name}.tmp
        rm {f2_name}.tmp
        rm {t}.tmp
        rm {s}.tmp
        rm {o1_name}.tmp
        rm {o2_name}.tmp
    done
fi

"""

    def __init__(self, qsub_option, script_dir):
        super(Bam2Fastq, self).__init__(qsub_option, script_dir)


