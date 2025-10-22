# -*- coding: utf-8 -*-
# @Time    : 2022/11/08 13:28
# @Author  : zzk
# @File    : demult.py
# Description:质粒根据yaml来分流。（新）
import glob
import gzip
import os.path
import time

import yaml

R_yaml = f'D:/zzk/NPdata/N110/BC13-N110/N110YZXI.yaml'
D_yaml = f'D:/zzk/NPdata/N110/BC13-N110/N110YZXD.yaml'
def r_input(yaml_file):
    # Check file format
    if not yaml_file.endswith('.yaml'):
        return
    with open(yaml_file, "r", encoding="utf8") as f:
        input_info = yaml.load(f, Loader=yaml.FullLoader)
    return input_info


def input_info(ref_yaml=R_yaml, demulti_yaml=D_yaml):
    """Extraction of excel information"""
    global ref_info, demul_info
    if ref_yaml:
        ref_info = r_input(ref_yaml)
    if demulti_yaml:
        demul_info = r_input(demulti_yaml)
    return ref_info, demul_info


def re_com(seq):
    return seq.replace("A", "B").replace("G", "D").replace("T", "A").replace("C", "G").replace("B",
                                                                                               "T").replace(
        "D", "C")[::-1]


def mkgreps():
    ref, dem = input_info(ref_yaml=D_yaml, demulti_yaml=D_yaml)
    for key in dem.keys():
        L_45 = dem[f"{key}"]["left_seqs"][:45]
        R_45 = dem[f"{key}"]["right_seqs"][-45:]
        grepss_L = [L_45[0 + 7 * i:17 + 7 * i] for i in range(round(45 / 7))]
        greps_L = [seq for seq in grepss_L if len(seq) == 17]
        gre_lwrite = ''
        for lw in greps_L:
            gre_lwrite += (lw + "\n")
        with open(f'D:/zzk/NPdata/N110/BC13-N110/greps/{key}-L.txt', 'w') as fl:
            fl.write(gre_lwrite)
            fl.close()

        grepss_R = [R_45[0 + 7 * i:17 + 7 * i] for i in range(round(45 / 7))]
        greps_R = [seq for seq in grepss_R if len(seq) == 17]
        gre_rwrite = ''
        for rw in greps_R:
            gre_rwrite += (rw + "\n")
        with open(f'D:/zzk/NPdata/N110/BC13-N110/greps/{key}-R.txt', 'w') as fr:
            fr.write(gre_rwrite)
            fr.close()


def read_fastq(fastq):
    ref, dem = input_info(ref_yaml=R_yaml, demulti_yaml=D_yaml)
    fq = gzip.open(fastq, 'rb')
    id_ = []
    with fq as f:
        a = 0
        while True:
            count_re = 0
            l1 = f.readline()
            if not l1:
                break
            # elif l1.decode('utf-8') == '\n':
            #     continue
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            read = l1 + l2 + l3 + l4
            # print(read)
            if isinstance(read, bytes):
                read = str(read, encoding='utf-8')
            read = tuple(filter(None, read.split('\n')))  # read[0]read[1]read[2]read[3]
            seq = read[1]  # 原序列
            # seq_q = read[3]
            seq_re = read[1].replace("A", "B").replace("G", "D").replace("T", "A").replace("C", "G").replace("B",
                                                                                                             "T").replace(
                "D", "C")[::-1]  # 反向互补序列
            # seq_re_q = seq_q[::-1]
            flag_repeat = []  # 检查交叉污染
            flag_repeat_num = 0
            for key in dem.keys():
                b, c, d = 0, 0, 0
                prior_L = b
                prior_R = c
                prior_LR = d

                L_45 = dem[f"{key}"]["left_seqs"][:45]
                R_45 = dem[f"{key}"]["right_seqs"][-45:]
                name_p = key
                grepss_L = [L_45[0 + 7 * i:17 + 7 * i] for i in range(round(45 / 7))]
                greps_L = [seq for seq in grepss_L if len(seq) == 17]

                grepss_R = [R_45[0 + 7 * i:17 + 7 * i] for i in range(round(45 / 7))]
                greps_R = [seq for seq in grepss_R if len(seq) == 17]
                # left
                flag_L = False
                seq_write, readQ = '', ''  # 存seq和对应的q值编码
                flagL_lis = [False, False, False, False]
                for grep in greps_L:
                    if grep in seq[:45]:
                        flag_L = True
                        flagL_lis[0] = True
                        seq_write = re_com(seq=seq)
                        readQ = read[3][::-1]
                    elif grep in seq[::-1][:45]:
                        count_re += 1
                        # print("re_L",key,count_re)
                        flag_L = True
                        flagL_lis[1] = True
                        seq_write = re_com(seq[::-1])
                        readQ = read[3]
                    elif grep in seq_re[:45]:
                        # print("1111111")
                        flag_L = True
                        flagL_lis[2] = True
                        seq_write = re_com(seq_re)
                        readQ = read[3]
                    elif grep in seq_re[::-1][:45]:
                        # print("re2_L",key, count_re)
                        flag_L = True
                        flagL_lis[3] = True
                        seq_write = re_com(seq_re[::-1])
                        readQ = read[3][::-1]
                    if flag_L == True:
                        break
                if flag_L:
                    seq_write = re_com(seq_write)
                    with open(f"D:/zzk/NPdata/N110/BC13-N110/Demultiplexed/{name_p}-L.fastq", 'a') as ff:
                        ff.write(read[0] + '\n' + seq_write + '\n' + read[2] + '\n' + readQ[::-1] + '\n')
                    b += 1
                seq_write, readQ = '', ''
                # Right
                flag_R = False
                flagR_lis = [False, False, False, False]
                for grep in greps_R:
                    if grep in seq[-45:]:
                        flag_R,flagR_lis[0] = True,True
                        seq_write,readQ = seq,read[3]
                    elif grep in seq[::-1][-45:]:
                        flag_R,flagR_lis[1] = True,True
                        seq_write,readQ = seq[::-1],read[3][::-1]
                    elif grep in seq_re[-45:]:
                        # print("re_R", key, count_re)
                        flag_R = True
                        flagR_lis[2] = True
                        seq_write = seq_re
                        readQ = read[3][::-1]
                    elif grep in seq_re[::-1][:45]:
                        # print("re2_R", key, count_re)
                        flag_R = True
                        flagR_lis[3] = True
                        seq_write = seq_re[::-1]
                        readQ = read[3]
                    if flag_R == True:
                        break
                if flag_R:
                    with open(f"D:/zzk/NPdata/N110/BC13-N110/Demultiplexed/{name_p}-R.fastq", 'a') as ff:
                        ff.write(read[0] + '\n' + seq_write + '\n' + read[2] + '\n' + readQ + '\n')
                    c += 1
                for ii in range(len(flagL_lis)):
                    if flagL_lis[ii] and flagR_lis[ii]:
                        with open(f"D:/zzk/NPdata/N110/BC13-N110/Demultiplexed/{name_p}-LR.fastq", 'a') as ff:
                            ff.write(read[0] + '\n' + seq_write + '\n' + read[2] + '\n' + readQ + '\n')
                        d += 1
                if prior_L - b > 0 or prior_R - c > 0 or prior_LR - d > 0:
                    flag_repeat.append(name_p)
                    flag_repeat_num += 1
            if len(flag_repeat) > 1:
                print("这里有交叉污染风险，同一个数据存在于多个fastq：", flag_repeat)
                id_.append(read[0])
            a += 1
            if a % 100000 == 0:
                print("已经处理了", a, "条数据,用时", round((time.time() - start) / 60, 2), "minutes.")
    return id_


def stats_fastq():
    file_fq = glob.glob(f"D:/zzk/NPdata/N110/BC13-N110/Demultiplexed/*.fastq")
    for file in file_fq:
        with open(file, 'r') as f:
            s = f.readlines()
            num = int(len(s) / 4)
            f.close()
        line = f"{file}抓取出来 {num} 条。"
        with open(f"D:/zzk/NPdata/N110/BC13-N110/Demultiplexed/NEW-Demultiplex_stats.txt", "a") as at:
            at.write(line + '\n')


if __name__ == '__main__':
    fastq = f'D:/zzk/NPdata/N110/BC13-N110/N110chop.fastq.gz'
    start = time.time()

    if not os.path.exists(f"D:/zzk/NPdata/N110/BC13-N110/greps"):
        os.mkdir(f"D:/zzk/NPdata/N110/BC13-N110/greps")
    # mkgreps()
    # id = read_fastq(fastq)
    stats_fastq()
    print(id)
    print("用时:", round((time.time() - start) / 60, 2))
