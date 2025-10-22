# -*- coding: utf-8 -*-
# @Time    : 2022/09/29 9:59
# @Author  : zzk
# @File    : fastq_visual.py
# Description:可视化。自动获取当前目录的fasta文件名，去找到对应的name.mmi、fastq,不过fastq名字要改改，每次用的不一样。

import glob
import os


def run(file):
    name_p = file[:-6]
    #cmd1 = f'minimap2 -x map-ont -d "{name_p}.mmi" "{file}"'
    cmd1 = f'bwa index "{file}"'
    os.system(cmd1)
    #cmd20 = f'minimap2 -ax map-ont "{name_p}.mmi" "{name_p}.fastq" > "{name_p}.sam"'
    cmd20 = f'bwa mem "{file}" "{name_p}.fastq" > "{name_p}.sam"'
    cmd210 = f'samtools view -bS "{name_p}.sam" > "{name_p}.bam"'
    cmd220 = f'samtools sort -O bam -o "{name_p}.sorted.bam" -T temp "{name_p}.bam"'
    cmd230 = f'samtools index "{name_p}.sorted.bam"'
    os.system(cmd20)
    os.system(cmd210)
    os.system(cmd220)
    os.system(cmd230)

