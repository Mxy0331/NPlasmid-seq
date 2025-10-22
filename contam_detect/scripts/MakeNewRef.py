# -*- coding: utf-8 -*-
# @Time    : 2022/12/20 14:03
# @Author  : zzk
# @File    : MakeNewRef.py
# Description:制作新的参考序列，主要针对sgRNA
import shutil


def run(file_fa, basic_20):
    newName = file_fa[:-6]
    shutil.copyfile(file_fa, f'{newName}-ori.fasta')
    content, head = "", ""
    with open(file_fa, 'r') as fa:
        s1 = fa.readline()
        head = s1.replace("\n", "")
        while s1:
            s1 = fa.readline()
            content += s1
        fa.close()
    content = content.replace("\n", "").replace(" ", "")  # 参考序列
    new_fa = ""
    n = 0
    for i in (content[:content.rfind(basic_20) + len(basic_20) + 10]):
        n += 1
        if n % 60 == 0:
            new_fa += f"{i}\n"
        else:
            new_fa += i

    # print(head+"\n"+new_fa)
    with open(file_fa,'w') as nfa:
        nfa.write(head+"\n"+new_fa)
        nfa.close()
    return file_fa,content[content.rfind(basic_20) -40:content.rfind(basic_20)-10]
# run("./data/S1263-2.fasta", "ATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGG")
