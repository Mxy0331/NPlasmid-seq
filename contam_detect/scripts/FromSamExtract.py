# -*- coding: utf-8 -*-
# @Time    : 2022/09/12 10:30
# @Author  : zzk
# @File    : FromSamExtract.py
# Description:读取Sam文件，解读CIGAR字符串。
# 普通质粒：根据CIGAR字符串，解析每个测序read跟参考序列的比对情况，将完整的数据放到normaldata.fastq，预警数据放到waring.fastq
# sgRNA：根据骨架定位独特碱基位置，取出独特碱基与参考比对，分normal与waring

import os
import nw as NW

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# start = time.time()

def read_cigar_(cigar_this_id_):
    num_lis = []
    flag_lis = []
    num = ''
    for j in range(len(cigar_this_id_)):
        if cigar_this_id_[j].isalpha():
            flag_lis.append(cigar_this_id_[j])
            num_lis.append(num)
            num = ''
        else:
            num += cigar_this_id_[j]
    len_re = 0
    len_read = 0
    for k in range(len(flag_lis)):
        if flag_lis[k] == 'M' or flag_lis[k] == 'D':
            len_re += int(num_lis[k])
        if flag_lis[k] == 'M' or flag_lis[k] == 'I':
            len_read += int(num_lis[k])
    return len_re, len_read


def read_nw_result(clustal_result):
    count_equal = 0  # 计算比例，求的是对上的总数
    seq_20_nt = clustal_result[1].replace("*", "N")
    count_20_ = 0
    count_continuous = 0  # 记录连续对上几个
    count_continuous_max = 0  # 记录连续对上的最大值
    seq_20_start, seq_20_end = 0, 0
    for cl in range(len(seq_20_nt)):
        if seq_20_nt[cl] != "N":
            count_20_ += 1
            if count_20_ == 1:
                seq_20_start = cl
            if count_20_ == 20:
                seq_20_end = cl
                break
    equal_list = []  # 记录所有连续比对个数
    for cll in range(seq_20_start, seq_20_end + 1):  # 只计算那20独特碱基的比对，不看sgRNA公共部分
        if clustal_result[0][cll] == clustal_result[1][cll]:
            count_continuous += 1

            if count_continuous_max < count_continuous:
                count_continuous_max = count_continuous
            count_equal += 1
        else:
            equal_list.append(count_continuous)
            count_continuous = 0
    equal_list.append(count_continuous)
    return seq_20_start, seq_20_end, equal_list, count_continuous_max, count_equal


def read_fasta(file):
    with open(file, 'r') as fa:
        con_fa = fa.readline()
        cont = ''
        while con_fa:
            if con_fa.startswith(">"):
                con_fa = fa.readline()
                continue
            else:
                cont += con_fa.replace("\n", "").replace(" ", "")
                con_fa = fa.readline()
    return cont.upper()


def find_seqs_by_id(reads_id, file, str):  # 存放无@的id
    name_p = file.split('\\')[-1][:-6]
    print(file)
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
            # else:
            #     with open(f"../929_result/{name_p}-normal.fastq", 'a') as f2:
            #         f2.write(con1 + con2 + con3 + con4)
            con1 = f1.readline()
            con2 = f1.readline()
            con3 = f1.readline()
            con4 = f1.readline()
    f1.close()


def find_seqs_by_20pos(infor, file, str):  # 存放无@的id
    name_p = file.split('\\')[-1][:-6]
    with open(file, 'r') as f1:
        con1 = f1.readline()
        con2 = f1.readline()
        con3 = f1.readline()
        con4 = f1.readline()
        while con1:
            id = con1.replace("\n", "").split(" ")[0].strip("@")
            if id in infor.keys():
                pos = infor[f'{id}']
                seq = con2[pos[0]:pos[1]] + "\n"
                q_seq = con4[pos[0]:pos[1]] + "\n"
                start_id = "@" + id + "\n"
                with open(f"./{name_p}-{str}.fastq", 'a') as f2:
                    f2.write(start_id + seq + con3 + q_seq)
                    f2.close()
            # else:
            #     with open(f"../929_result/{name_p}-normal.fastq", 'a') as f2:
            #         f2.write(con1 + con2 + con3 + con4)
            con1 = f1.readline()
            con2 = f1.readline()
            con3 = f1.readline()
            con4 = f1.readline()
    f1.close()


def read_fastq(file_fq):
    with open(file_fq, 'r') as fq:
        cont_fastq_ = fq.readlines()
        cont_fastq = []
        total_re = int(len(cont_fastq_) / 4)
        for o in range(len(cont_fastq_)):
            cont_fastq.append(cont_fastq_[o].replace("\n", "").split(" ")[0])
        fq.close()
    return cont_fastq


# 首先要知道reference的长度和read的长度。假设都知道了read_len_real,reference_len_real
def run(file_fa, file_fq, file_sam, seq_20nt):  # file_fa, file_fq, file_sam
    reference_len_real = len(read_fasta(file_fa))
    cont_fastq = read_fastq(file_fq=file_fq)
    record_20_pos = {}
    count_waring = 0
    count_waring_id,normal_id = [],[]
    count_abnormal, list_count_abnormal = 0, []
    # sam文件中的CIGAR要着重理解
    with open(f'{file_sam}') as f:  # 解析sam文件
        '''
        将连续同一个read的比对信息的cigar接到一起，去掉其中的S和H，计算长度？
        '''
        cont = f.readlines()  # 将sam文件所有内容拿出来
        id_old = ""  #
        i, count_all = 0, 0
        while i < len(cont) and cont[i]:  # 遍历sam文件的内容
            row = cont[i].split('\t')
            id, loc, cigar = row[0], row[1], row[3]
            the_cigar = {}  # 存放某id对应的所有比对起始位置和cigar
            read_real_seq = cont_fastq[cont_fastq.index(f"@{id}") + 1]  # 取出测序seq
            read_len_real = len(cont_fastq[cont_fastq.index(f"@{id}") + 1])  # 遍历的这个read的长度

            the_same_id_cigar, the_same_id_loc = [], []
            if id != id_old:  # 如果这个id跟前一个id不一样
                the_cigar[int(f'{loc}')] = cigar
            else:
                break
            for k in range(i + 1, len(cont)):  # 从i+1开始遍历后边的所有内容，找到与i是同一条read的内容
                # 每遍历一个就取到所有的关键信息
                row_same = cont[k].split('\t')
                id_same, loc_same, cigar_same = row_same[0], row_same[1], row_same[3]
                if id_same != id:  # 如果连续两个id不一样，跳过这里
                    break
                else:  # 同一个read的比对信息，就把pos和cigar添加到the_same_id_cigar和the_same_id_loc
                    id_old = id_same
                    # loc_old = loc_same
                    the_cigar[int(f'{loc_same}')] = cigar_same
            # 比对信息  the_cigar = {“1”：cigar1，“1839”：cigar2}

            # -----------start-----------------------
            blast_re_len = 0  # 参考序列比对区间
            blast_read_len = 0  # reads比对区间
            M_num, I_num, D_num = 0, 0, 0
            # cigar_length = 0
            # cigar_flag_length = 0
            for key_ori, value_ori in the_cigar.items():  # key_ori是位置，value_ori是CIGAR
                numk = ''
                flag_lis_key = []  # 记录MDISH
                num_lis_key = []  # 记录MDISH对应的num

                for jk in range(len(value_ori.replace("\n", ""))):  # 遍历这个CIGAR
                    if value_ori[jk].isalpha():
                        flag_lis_key.append(value_ori[jk])
                        num_lis_key.append(int(numk))
                        numk = ''
                    else:
                        numk += value_ori[jk]
                for kj in range(len(flag_lis_key)):
                    if flag_lis_key[kj] == "M":
                        blast_re_len += num_lis_key[kj]
                        blast_read_len += num_lis_key[kj]
                        M_num += num_lis_key[kj]
                    if flag_lis_key[kj] == "D":
                        blast_re_len += num_lis_key[kj]
                        D_num += num_lis_key[kj]
                    if flag_lis_key[kj] == "I":
                        blast_read_len += num_lis_key[kj]
                        I_num += num_lis_key[kj]
            re_result = round(blast_re_len / reference_len_real, 4)
            abnormal = False
            if re_result < 0.95 or re_result > 1.02:  # 对分段的也统计了比对长度百分比
                abnormal = True
            # -------------end----------------

            for key in sorted(the_cigar):  # 根据CIGAR对应的起始位置排序，小的在前。
                # print(key,type(key))
                the_same_id_cigar.append(the_cigar[key])
                the_same_id_loc.append(key)
            '''
            11.04 拿出每一条sgRNA read的那20bp
            '''
            s_start = 0  # CIGAR开头S对应的数值
            blast_re_, blast_read_ = 0, 0
            if not flag_lis_key:
                i += 1
                continue
            if flag_lis_key[0] == "S":
                s_start = num_lis_key[0]

            for kj in range(len(flag_lis_key)):
                if flag_lis_key[kj] == "M":
                    blast_re_ += num_lis_key[kj]
                    blast_read_ += num_lis_key[kj]
                if flag_lis_key[kj] == "D":
                    blast_re_ += num_lis_key[kj]
                if flag_lis_key[kj] == "I":
                    blast_read_ += num_lis_key[kj]
            the_rest_relen = reference_len_real - blast_re_ - int(loc) + 1

            if blast_read_ > 0 or blast_re_ > 0:  #
                start_20 = blast_read_ + s_start + the_rest_relen - 34
                end_20 = blast_read_ + s_start + the_rest_relen
                record_20_pos[f"{id}"] = [start_20, end_20]
                refer = f'NNCG{seq_20nt}GTTTCAGAGNNNNNNNNNN'
                clustal_result = NW.dynamic_edit(read_real_seq[start_20:end_20][3:],
                                                 refer[:len(read_real_seq[start_20:end_20])][3:])
                seq_20_start, seq_20_end, equal_list, count_continuous_max, count_equal = read_nw_result(clustal_result)

                if read_len_real < 4200:  # 排除多倍长度，因为多倍长度发生率很少，污染载体基本不会产生多倍长度
                    # 其实也可以排除reference_len_real*110%这种异常
                    right_rate = round(count_equal / (seq_20_end - seq_20_start + 1), 3)
                    if right_rate < 0.6 and count_continuous_max < 6 and (
                            equal_list.count(4) + equal_list.count(5)) < 2:
                        count_waring += 1
                        count_waring_id.append(id)
                        # print("预警！！")
                        if abnormal:
                            count_abnormal += 1
                            list_count_abnormal.append(id)  # 异常
            # ------------End-------------
            if id not in count_waring_id:
                normal_id.append(id)
            count_all += 1
            if i == k:
                break
            else:
                i = k
        logger.info(f"{file_fq}: after minimap,the reads are {count_all}")
    f.close()

    find_seqs_by_20pos(record_20_pos, f"{file_fq}", "20")
    find_seqs_by_id(count_waring_id, f"{file_fq[:-6]}-20.fastq", "Waring")
    find_seqs_by_id(normal_id, f"{file_fq[:-6]}.fastq", "NormalData")
    os.remove(f"{file_fq[:-6]}-20.fastq")
    logger.info(f"waring reads are {len(count_waring_id)}")
    logger.info(f"waring reads has saved in '{file_fq[:-6]}-20-waring.fastq'")  # 对这个文件聚类
    return f'{file_fq[:-6]}-20-waring.fastq',count_all


def conmmonRun(refervalue, readvalue, file_fa, file_fq, file_sam):
    reference_len_real = len(read_fasta(file_fa))  # 参考序列长度
    cont_fastq = read_fastq(file_fq=file_fq)  # 测序数据
    # record_20_pos = {}
    count_waring = 0
    count_waring_id,normal_id = [],[]
    waring_blast_result = {}
    count_abnormal, count_po_1263, list_count_abnormal, list_count_po_1263 = 0, 0, [], []

    with open(f'{file_sam}', 'r') as f:  # 解析sam文件
        '''
        2022-12-23将连续同一个read的比对信息的cigar接到一起，去掉其中的S和H，计算长度
        正常minimap比对是不会重复比对的，上述可行
        '''
        cont = f.readlines()  # 将sam文件所有内容拿出来
        id_old = ""  #
        i, count_all = 0, 0
        while i < len(cont) and cont[i]:  # 遍历sam文件的内容
            row = cont[i].split('\t')
            id, loc, cigar = row[0], row[1], row[3]
            the_cigar = {}  # 存放某id对应的所有比对起始位置和cigar
            read_real_seq = cont_fastq[cont_fastq.index(f"@{id}") + 1]  # 取出测序seq
            read_len_real = len(cont_fastq[cont_fastq.index(f"@{id}") + 1])  # 遍历的这个read的长度

            the_same_id_cigar, the_same_id_loc = [], []
            if id != id_old:  # 如果这个id跟前一个id不一样
                the_cigar[int(f'{loc}')] = cigar
            else:
                break
            for k in range(i + 1, len(cont)):  # 从i+1开始遍历后边的所有内容，找到与i是同一条read的内容
                # 每遍历一个就取到所有的关键信息
                row_same = cont[k].split('\t')
                id_same, loc_same, cigar_same = row_same[0], row_same[1], row_same[3]
                if id_same != id:  # 如果连续两个id不一样，跳过这里
                    break
                else:  # 同一个read的比对信息，就把pos和cigar添加到the_same_id_cigar和the_same_id_loc
                    id_old = id_same
                    # loc_old = loc_same
                    the_cigar[int(f'{loc_same}')] = cigar_same
            # 比对信息  the_cigar = {“1”：cigar1，“1839”：cigar2}

            # -----------start-----------------------
            blast_re_len = 0  # 参考序列比对区间
            blast_read_len = 0  # reads比对区间
            M_num, I_num, D_num = 0, 0, 0

            for key_ori, value_ori in the_cigar.items():  # key_ori是位置，value_ori是CIGAR
                numk = ''
                flag_lis_key = []  # 记录MDISH
                num_lis_key = []  # 记录MDISH对应的num

                for jk in range(len(value_ori.replace("\n", ""))):  # 遍历这个CIGAR
                    if value_ori[jk].isalpha():
                        flag_lis_key.append(value_ori[jk])
                        num_lis_key.append(int(numk))
                        numk = ''
                    else:
                        numk += value_ori[jk]
                for kj in range(len(flag_lis_key)):
                    if flag_lis_key[kj] == "M":
                        blast_re_len += num_lis_key[kj]
                        blast_read_len += num_lis_key[kj]
                        M_num += num_lis_key[kj]
                    if flag_lis_key[kj] == "D":
                        blast_re_len += num_lis_key[kj]
                        D_num += num_lis_key[kj]
                    if flag_lis_key[kj] == "I":
                        blast_read_len += num_lis_key[kj]
                        I_num += num_lis_key[kj]
            re_result = round(blast_re_len / reference_len_real, 4)
            read_result = round(blast_read_len / read_len_real, 4)
            if re_result < refervalue or read_result < readvalue:  # 比对区域小于某个数值的看作异常数据
                count_waring_id.append(id)
                waring_blast_result[f'{id}'] = [read_len_real,re_result,read_result]
            else:
                normal_id.append(id)

            # -------------end----------------
            count_all += 1
            if i == k:
                break
            else:
                i = k
        logger.info(f"{file_fq}: after minimap,the reads are {count_all}")
    f.close()

    # find_seqs_by_20pos(record_20_pos, f"{file_fq}", "20")
    find_seqs_by_id(count_waring_id, f"{file_fq[:-6]}.fastq", "Waring")
    find_seqs_by_id(normal_id, f"{file_fq[:-6]}.fastq", "NormalData")
    logger.info(f"waring reads are {len(count_waring_id)}")
    logger.info(f"waring reads has saved in '{file_fq[:-6]}-Waring.fastq'")  # 对这个文件聚类

    return f'{file_fq[:-6]}-Waring.fastq',waring_blast_result,count_all


def main(file_fa, file_fq, file_sam, seq_20nt=None, readvalue=None, refervalue=None):
    if seq_20nt:
        cluster_fq,allData = run(file_fa=file_fa, file_fq=file_fq, file_sam=file_sam, seq_20nt=seq_20nt)
        return cluster_fq,allData
    else:
        cluster_fq,waring_blast_result,allData = conmmonRun(file_fa=file_fa, file_fq=file_fq, file_sam=file_sam, readvalue=readvalue,
                                refervalue=refervalue)
        return cluster_fq,waring_blast_result,allData
