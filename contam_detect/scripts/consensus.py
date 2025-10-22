# -*- coding: utf-8 -*-
# @Time    : 2022/12/19 14:39
# @Author  : zzk
# @File    : consensus.py
# Description:提共识序列
import glob
import os
import time

# cmd0 = f'conda init'
# os.system(cmd0)
'''
把consen中的文件格式化，然后提共识序列
'''


def consensus():
    '''
    把consen下的所有fastq名字格式化，然后提共识序列
    '''
    file_remate = glob.glob("./consen/*.fastq")
    sub = '.fastq'
    i = 0
    for file in file_remate:
        file_name = file.split("/")[-1].strip(sub)
        file_rename = f'./consen/umi{i}_{file_name}_bins{sub}'
        os.rename(file, file_rename)
        i += 1
    time.sleep(2)
    cmd2 = f"longread_umi consensus_racon -d ./consen -o ./racon -r 3 -t 16 -p 'map-ont'"
    os.system(cmd2)
    return "./racon/consensus_racon.fa"


def run(left_seq, file_fa):
    file_con = consensus()
    greps_left = [left_seq[-15:], left_seq[-18:-3], left_seq[-21:-6], left_seq[-24:-9], left_seq[-27:-12]]
    print(greps_left)

    dic = {}  # 添加一个dic[f'reference-{id}']=参考sgRNA
    f = open(file_fa)
    faName = f.readline().replace("\n", "")
    ReferenceSeq = ''.join(f.readlines()).replace('\n', '')
    with open(file_con, 'r') as fc:
        id = fc.readline().replace("\n", "").split("ubs=")[0]
        seq = fc.readline().replace("\n", "")
        while id:
            nu, seq_part = 0, ''
            for seql in greps_left:
                nu += 1
                if seq.rfind(seql) != -1:
                    start, seq_part, re_start = seq.rfind(seql) + 15, seql, ReferenceSeq.rfind(seql)
                    break
            print(start)
            dic[f'{id}'] = seq[start:start + nu * 15 + 40]
            dic[f'reference-{id}'] = ReferenceSeq[re_start:re_start + nu * 15 + 40]
            id = fc.readline().replace("\n", "").split("ubs=")[0]
            seq = fc.readline().replace("\n", "")
    return dic
