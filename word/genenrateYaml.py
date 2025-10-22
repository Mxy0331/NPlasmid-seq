# -*- coding: utf-8 -*-
# @Time    : 2022/07/14 15:16
# @Author  : zzk
# @File    : generateYaml.py
# Description: 9.15 N104修改：增加了对P2151A P2151B的判断
#2023-3.13
#2023-3-21 增加L R单侧   -zzk
#2023-3-28   word取begin\end改

import sys
import time
import pandas as pd
import yaml
import os
import xlrd
import glob
import docx


def find_loca(res, name_p):
    temp_location = 0
    for line in res:
        # print(line)
        info = line.split('	')
        print(info)
        if info[0] == name_p:
            f = info[3]
            if f == '+':
                temp_location = int(info[5]) - 2
            elif f == '-':
                temp_location = int(info[4]) + 3
            print(temp_location)
    return temp_location

def rev_complement(seq):
    """Get reverse complement of a DNA sequence"""
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def generate_yaml(name_p, RefInfo_name, DemultiInfo_name, seq_, strand, sgRNA, uniq_seq):
    '''
    需要copy过来的TempPl.txt,放在当前目录
    需要输入name_P: 当前处理的序列号
    需要输入RefInfo_name: 参考序列yaml文件名
    需要输入DemultiInfo_name: 分流yaml文件名
    需要输入seq_: “seqkit locate TempP1569.fa -p CATGAAGCGCCACCTTTG” 中的这个序列
    '''
    # name_p = f'P1677'
    # RefInfo_name = f'N98RefInfo'
    # DemultiInfo_name = f'N98DemultiInfor'
    # seq_ = f'CATGAAGCGCCACCTTTG'
    seq_ = seq_.replace("\n", "")#creSeq
    cmd_to_fa = f'seqkit fx2tab Temp{name_p}Pl.txt  | seqkit tab2fx | seqkit seq -u | seqkit replace -p " |-|\d" -s  > Temp{name_p}.fa'
    os.system(cmd_to_fa)
    # with open(f'Temp{name_p}.fa','r') as f1:
    #     a = f1.readlines()
    #     s=''
    #     for lin in range(len(a)):
    #         if a!=0:
    #             s += a[lin].replace("\n","")
    #     # print("那个东西：\n",s)
    cmd_check_num = f'seqkit locate Temp{name_p}.fa -p {seq_}'
    res = os.popen(cmd_check_num)
    temp_loc = find_loca(res, name_p)  # 得到切点。
    cmd_generate_fa = f'seqkit restart Temp{name_p}.fa -i {temp_loc} > {name_p}fa.txt'

    os.system(cmd_generate_fa)

    '''
    根据之前的样子定义了yaml文件的一些东西。
    '''
    rev_sgRNA = ''
    if strand == '-':
        rev_sgRNA = rev_complement(sgRNA)

    dic_sample_LR = {
        name_p: {"reference_id": name_p, "sgRNA": sgRNA, "rev_com_sgRNAseq": rev_sgRNA, "left_seqs": '', "right_seqs": '',
                 "BCprimer_F": '', "BClen_F": '',
                 "BCprimer_R": '', "BClen_R": '', "unique_sequence": '', "strand": strand}}
    dic_sample_L = {
        f'{name_p}-L': {"reference_id": f"{name_p}", "sgRNA": sgRNA, "rev_com_sgRNAseq": rev_sgRNA, "left_seqs": '',
                        "right_seqs": '', "BCprimer_F": '', "BClen_F": '',
                        "BCprimer_R": '', "BClen_R": '', "unique_sequence": '', "strand": strand}}
    dic_sample_R = {
        f'{name_p}-R': {"reference_id": f"{name_p}", "sgRNA": sgRNA, "rev_com_sgRNAseq": rev_sgRNA, "left_seqs": '',
                        "right_seqs": '', "BCprimer_F": '', "BClen_F": '',
                        "BCprimer_R": '', "BClen_R": '', "unique_sequence": '', "strand": strand}}
    dic_sample_ref = {name_p: {"sequence": ''}}
    with open(f'{name_p}fa.txt', 'r', encoding='utf-8') as ref:
        seq_ref = ref.readlines()
        seq_ref_all = ""
        for i in range(len(seq_ref)):
            seq_ref_all += seq_ref[i]
        seq_ref_all = seq_ref_all.replace(">", "").replace(name_p, "").replace("\t", "").replace("\n", "").replace(" ",
                                                                                                                   "")
        print("reference length:", len(seq_ref_all))
        # print(seq_ref_all)
        dic_sample_ref[name_p]['sequence'] = seq_ref_all  # N98RefInfo.yaml的dict就好了
    with open(f'{RefInfo_name}.yaml', 'a', encoding='utf-8') as f:
        yaml.dump(dic_sample_ref, f, default_style=False, encoding='utf-8', allow_unicode=True)

    dic_sample_LR[name_p]["left_seqs"] = seq_ref_all[:40]
    dic_sample_LR[name_p]["right_seqs"] = seq_ref_all[-40:]
    dic_sample_L[f'{name_p}-L']["left_seqs"] = seq_ref_all[:40]
    dic_sample_R[f'{name_p}-R']["right_seqs"] = seq_ref_all[-40:]
    # dic_sample_LR[name_p]["unique_sequence"] = uniq_seq

    with open(f'{DemultiInfo_name}.yaml', 'a', encoding='utf-8') as f:
        yaml.dump(dic_sample_LR, f, default_style=False, encoding='utf-8', allow_unicode=True)
    with open(f'{DemultiInfo_name}.yaml', 'a', encoding='utf-8') as f:
        yaml.dump(dic_sample_L, f, default_style=False, encoding='utf-8', allow_unicode=True)
    with open(f'{DemultiInfo_name}.yaml', 'a', encoding='utf-8') as f:
        yaml.dump(dic_sample_R, f, default_style=False, encoding='utf-8', allow_unicode=True)

def read_docx(file_doc, name_p):
    '''
    word都是以序列号结尾
    读取docx，生成TempPl.txt
    调用这个返回调整好的序列，直接可以放到Temptxt中
    '''
    # 获取文档对象
    file = docx.Document(file_doc)
    seq = ''
    doc_ = []
    count, startofseq, endofseq = 0, 0, 0

    doc_ = []
    count, startofseq, endofseq = 0, 0, 0
    found_begin = False

    for para in file.paragraphs:
        count += 1
        if 'startofseq' in para.text.lower():
            startofseq = count
            found_begin = True
        if found_begin and 'endofseq' in para.text.lower():
            print(para.text.lower())
            print(f"==={count}")
            endofseq = count
            break
        doc_.append(para.text)

    for i in range(startofseq, endofseq - 1):  # .isalpha()
        for j in range(len(doc_[i])):
            if doc_[i][j].isalpha():
                seq += doc_[i][j]
    seq_plasmid = seq.replace("\n", "").replace(" ", "").upper().replace("+", "")
    len_mid = int(len(seq_plasmid) / 2)
    uniq_seq = seq_plasmid[len_mid:len_mid + 201]
    with open(f'Temp{name_p}Pl.txt', 'a') as f:
        f.write(f">{name_p}\n")
        f.write(seq_plasmid + '\n')

    return uniq_seq, seq_plasmid


if __name__ == '__main__':

    print("当前目录只能包含一个有“Pre”的xlsx文件！！！！")
    all_file_doc = glob.glob("*.docx")  # 取到所有的word文件
    if len(all_file_doc) == 0:
        print("当前目录没有word文件，程序三秒后关闭！")
        time.sleep(3)
        sys.exit()

    file_list_xlsx = glob.glob("." + "/*.xlsx")
    file_pre_name = ""
    for file in file_list_xlsx:
        if "质粒-sgRNA" in file:
            file_pre_name = file
    df = pd.read_excel(file_pre_name, engine='openpyxl')  # 使用pandas读取excel文件
    rows_count = len(df)

    print("如果多次运行处理同一批数据，需要将之前生成的yaml删除，原因是数据时依次不覆盖追加！")
    RefInfo_name = 'refinfo'
    DemultiInfo_name = 'Demuinfo'

    for i in range(0, rows_count):
        # 遍历要处理的质粒
        row_data = df.loc[i].values.tolist()  # 使用pandas获取行数据
        print(row_data[0], row_data[1])  # 0是质粒序号   1是切割位点序列
        uniq_seq = ''
        name_p = row_data[0].replace("\n", "").replace(" ", "") #改为自动识别名字和序列
        for file_doc in all_file_doc:
            name_p_no_AB = name_p.strip("A").strip('B')#只要名字和参考序列部分匹配
            if file_doc.endswith(f'{name_p}.docx') or file_doc.startswith(f'{name_p}'):
                print("name_p:",name_p)
                print("doc:",file_doc)
                uniq_seq, seq_plasmid = read_docx(file_doc=file_doc, name_p=name_p)
                break
        seq_ = row_data[3].replace("\n", "").replace(" ", "")
        sgRNA = row_data[3].replace("\n", "").replace(" ", "")
        strand = row_data[5].replace("\n", "").replace(" ", "")
        generate_yaml(name_p=name_p, RefInfo_name=RefInfo_name, DemultiInfo_name=DemultiInfo_name, seq_=seq_, strand=strand, sgRNA=sgRNA,
                      uniq_seq=uniq_seq)  # 生成了yaml

    all_file_t = os.listdir('.')
    print(all_file_t)
    for file in all_file_t:
        if file.endswith(".txt") or file.endswith(".fa"):
            os.remove(file)
    print("Finished")

    # yaml_or_demu = input("是否根据生成的yaml对数据进行分流？（y or n）：\n")
    # if yaml_or_demu == 'y':
    #     print("注：当前目录下有greporeseq代码、去接头后的数据")
    #     print("***************************************")
    #     chop_name = input("请输入去接头后的文件名（.fastq.gz之前的名）：\n")
    #     cmd_demu = f'python greporeseq1_4/greporeseq.py demultiplex -n {chop_name}.fastq.gz -d {DemultiInfo_name}.yaml'
    #     os.system(cmd_demu)

    # gen_reference_fasta = input("是否生成reference（y or n）：\n")
    # if gen_reference_fasta == 'y':
    #     print("注：当前目录下需有info.yaml")
    #     print("***************************************")
    #     cmd_info = f'python greporeseq1_4/greporeseq.py mkreference -r {RefInfo_name}.yaml'
    #     os.system(cmd_info)
    # print("over............")
