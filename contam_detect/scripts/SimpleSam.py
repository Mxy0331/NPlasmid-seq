# -*- coding: utf-8 -*-
# @Time    : 2022/09/12 10:30
# @Author  : zzk
# @File    : read_sam.py
# Description:简化（过滤）sam文件内容
import os

def run(file):
    if os.path.exists(f'{file[:-4]}-simple.sam'):
        return f'{file[:-4]}-simple.sam'
    with open(f'{file}') as f:
        cont = f.readline()
        while cont:
            if cont.startswith("@"):
                cont = f.readline()
                continue
            cont_lis = cont.split("\t")
            with open(f'{file[:-4]}-simple.sam', 'a') as ci:
                ci.write(cont_lis[0] + '\t' + cont_lis[3] + '\t' + cont_lis[4] + '\t' + cont_lis[5] + '\n')
                ci.close()
            cont = f.readline()
        f.close()
    os.remove(f'{file}')
    return f'{file[:-4]}-simple.sam'
