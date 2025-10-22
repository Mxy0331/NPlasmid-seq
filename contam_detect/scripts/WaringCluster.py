# -*- coding: utf-8 -*-
# @Time    : 2022/12/01 9:15
# @Author  : zzk
# @File    : waring-cluster.py
# Description:
import os.path
import shutil

import nw as NW
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def find_seqs_by_id(reads_id, file, str):  # 存放无@的id
    name_p = file.split('\\')[-1][:-6]
    with open(file, 'r') as f1:
        con1 = f1.readline()
        con2 = f1.readline()
        con3 = f1.readline()
        con4 = f1.readline()
        while con1:
            id = con1.replace("\n", "").split(" ")[0].strip("@")
            if id in reads_id:
                with open(f"./{name_p}-{str}.fastq", 'a') as f2:
                    f2.write(con1 + con2 + con3 + con4)
                    f2.close()
            con1 = f1.readline()
            con2 = f1.readline()
            con3 = f1.readline()
            con4 = f1.readline()
    f1.close()


def read_waring_fq(file):
    dic_waring_fq = {}
    with open(file, 'r') as f:
        s1 = f.readline()
        s2 = f.readline()
        s3 = f.readline()
        s4 = f.readline()
        while s1:
            id = s1.replace("\n", "").split(" ")[0].strip("@")
            dic_waring_fq[f'{id}'] = s2.replace("\n", "")
            s1 = f.readline()
            s2 = f.readline()
            s3 = f.readline()
            s4 = f.readline()
    return dic_waring_fq


def commonRun(waring_fq, waring_blast_result, top_n):
    '''
    根据长度、比对区域相似性进行聚类。
    :param waring_fq: 测序预警data
    :param waring_blast_result: 存放  key-预警id,value-(测序长度、CIGAR参考比对范围和测序数据比对范围 )
    :return: 聚类结果
    '''
    print("22222", waring_fq)
    cluster_name = []
    name = waring_fq[:-13]
    dic_classes = {}  # 存放聚的类，key为0，1，2，3，，，，value为[{id1:dic1}，{id2:dic2}......]
    num_process = 0
    class_num = 1  # 类别数量
    for key, value in waring_blast_result.items():  # id=key
        Difference_min = 1  # 最小差异
        flag_class = -1  # 记录处理过程中value属于哪一类
        num_process += 1
        if num_process == 1:
            dic_classes[f'{0}'] = [{f'{key}': value}]
            continue
        for key_classes, value_classes in dic_classes.items():  # 遍历每一个类别 key_classes类别序号，value_classes-[{id1:dic1}，{id2:dic2}......]
            '''
            两两比较，在一定阈值范围内，并且是差异最小的归为一类，确保不重复
            '''
            for value_class in value_classes[0].values():  # 取当前类的第一个seq
                seqLenDeff = 1 - abs((value[0] - value_class[0]) / max(value[0], value_class[0]))  # reads长度差异
                # reBlastDeff = abs(value[1] - value_class[1])
                # readBlastDeff = abs(value[2] - value_class[2])
                Difference_ = abs(value[1] - value_class[1]) + abs(value[2] - value_class[2])
                if seqLenDeff > 0.95 and Difference_ < 0.05:
                    if Difference_min > Difference_:
                        flag_class = key_classes
                break
        if flag_class == -1:
            dic_classes[f'{class_num}'] = [{f'{key}': value}]
            class_num += 1
        else:
            dic_classes[f'{flag_class}'].append({f'{key}': value})
    # print(dic_classes)
    dic_classes_num = {}
    with open("cluster.txt", 'w') as re:
        for key_, value_ in dic_classes.items():
            re.write(f"The class {key_} has {len(value_)} reads.\n")
            dic_classes_num[f'{key_}'] = len(value_)
        re.close()
    dic_classes_num_sort = sorted(dic_classes_num.items(), key=lambda x: x[1], reverse=True)
    if len(dic_classes_num_sort) < top_n:
        top_n = len(dic_classes_num_sort)
    logger.info(f"the top {top_n}:")
    for tn in range(top_n):
        logger.info(f"第{dic_classes_num_sort[tn][0]}类共有{dic_classes_num_sort[tn][1]}reads")
        id = []
        for kk in range(len(dic_classes[f'{dic_classes_num_sort[tn][0]}'])):
            for key_id in dic_classes[f'{dic_classes_num_sort[tn][0]}'][kk].keys():
                id.append(key_id)
        print("11111", name)
        find_seqs_by_id(id, f'{name}.fastq', f'cluster-top-{tn}')
        logger.info(f"the reads has saved to '{name}-cluster-top-{tn}.fastq'")
        pa = f'{name}-cluster-top-{tn}.fastq'.split("/")[-1]
        shutil.copyfile(f'{name}-cluster-top-{tn}.fastq', f'./consen/{pa}')
        cluster_name.append(f'{name}-cluster-top-{tn}')
    return cluster_name


def run(waring_fq, top_n):
    waring_fq = waring_fq
    name = waring_fq[:-16]
    dic_waring_fq = read_waring_fq(waring_fq)  # 存放所有预警id及其序列

    dic_classes = {}  # 存放聚的类，key为0，1，2，3，，，，value为[{id1:seq1}，{id2:seq2}......]
    num_process = 0
    class_num = 1  # 类别数量
    for key, value in dic_waring_fq.items():  # 处理每一个seq，与每一类的第一个进行比较匹配度
        match_max = 0  # 最大匹配度
        flag_class = -1  # 记录处理过程中value属于哪一类
        num_process += 1

        if num_process == 1:
            dic_classes[f'{0}'] = [{f'{key}': value}]
            continue
        # print(key, value)
        for key_classes, value_classes in dic_classes.items():  # 遍历每一个类别
            # print(value_classes)
            for value_class in value_classes[0].values():  # 取当前类的第一个seq
                clustal_result = NW.dynamic_edit(value[:25], value_class[:25])  # 将此value与每个类别中的第一个进行比较（以第一个为代表序列）
                break  # clustal_result就是与该类的比对结果
            count_match, count_continuous, count_continuous_max = 0, 0, 0
            equal_list = []
            for i in range(len(clustal_result[1])):
                if clustal_result[0][i] == clustal_result[1][i]:
                    count_match += 1
                    count_continuous += 1
                    if count_continuous_max < count_continuous:
                        count_continuous_max = count_continuous
                else:
                    equal_list.append(count_continuous)
                    count_continuous = 0
            equal_list.append(count_continuous)
            match = round(count_match / len(clustal_result[1]), 2)

            if match > 0.6 and match > match_max and count_continuous_max >= 5:
                match_max = match
                flag_class = key_classes  # 该seq在当前处理过程中属于key_classes类，记录下类别

        # 将seq加入到match_max对应的类中，也就是flag_class代表的类,如果flag_class=-1，就新加一个类。
        if flag_class == -1:
            dic_classes[f'{class_num}'] = [{f'{key}': value}]
            class_num += 1
        else:
            dic_classes[f'{flag_class}'].append({f'{key}': value})

    dic_classes_num = {}
    with open("cluster.txt", 'w') as re:
        for key_, value_ in dic_classes.items():
            re.write(f"The class {key_} has {len(value_)} reads.\n")
            dic_classes_num[f'{key_}'] = len(value_)
        re.close()
    dic_classes_num_sort = sorted(dic_classes_num.items(), key=lambda x: x[1], reverse=True)
    if len(dic_classes_num_sort) < top_n:
        top_n = len(dic_classes_num_sort)
    logger.info(f"the top {top_n}:")
    for tn in range(top_n):
        logger.info(f"第{dic_classes_num_sort[tn][0]}类共有{dic_classes_num_sort[tn][1]}reads")
        id = []
        for kk in range(len(dic_classes[f'{dic_classes_num_sort[tn][0]}'])):
            for key_id in dic_classes[f'{dic_classes_num_sort[tn][0]}'][kk].keys():
                id.append(key_id)
        find_seqs_by_id(id, f'{name}.fastq', f'cluster-top-{tn}')
        logger.info(f"the reads has saved to '{name}-cluster-top-{tn}.fastq'")
        pa = f'{name}-cluster-top-{tn}.fastq'.split("/")[-1]
        shutil.copyfile(f'{name}-cluster-top-{tn}.fastq', f'./consen/{pa}')
