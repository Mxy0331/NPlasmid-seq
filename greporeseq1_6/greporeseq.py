# -*- coding: utf-8 -*-
# @Time         : 2021/7/27
# @Author       : Siang_Li
# @File         : greporeseq_ori.py
# @Software     : PyCharm
# @E-mail       : lisiang.darcy@gmail.com
# @Description  : a pipeline for analysing ONT sequencing data

import os  #提供了丰富的方法来处理文件和目录，同时提供了与操作系统进行交互的接口，如读取环境变量、操作文件和目录等。
import argparse #用于编写用户友好的命令行接口。它能处理程序参数，自动生成帮助和使用消息，并在用户给程序提供无效信息时发出错误。
import re #用于处理正则表达式的模块，提供了对正则表达式的支持。
import textwrap  #用于进行文本的换行操作。
import time  #提供了时间相关的功能，包括获取当前时间、操作时间和日期、实现延迟等。
import subprocess  #用于产生子进程，并连接到它们的输入/输出/错误管道，获取它们的返回码。
import glob  #用于在Unix样式的路径名模式扩展中生成文件列表。
from tkinter.messagebox import NO  #tkinter用于创建图形用户界面(GUI)。其中，messagebox.NO是一个常量，用于代表否定的消息框返回值；tix.Tree是一个类，用于创建和操作Tree控件。
from tkinter.tix import Tree
import random
from unicodedata import name #用于访问Unicode字符数据库。name函数可以返回指定Unicode字符的名称。

import log
import visualization as vis
import read_input as ri
from read_fastq import read_fastq
# 首先, 我们创建一个自定义的日志记录器，命名为 'root'
# 'createCustomLogger' 是一个可能在 'log' 模块中定义的函数，用于创建自定义日志记录器
# 其目的是为了记录和跟踪程序运行过程中的事件，提供调试信息，以及记录可能的错误和异常情况
logger = log.createCustomLogger('root')


# 定义一个名为 'rm_dup_none' 的函数，该函数接受一个列表 'l' 作为参数
# 函数首先通过 'set' 去除列表中的重复元素，然后使用 'sorted' 按原列表的顺序排序
# 使用 'filter' 函数删除所有的 'None' 元素
# 最后，函数返回处理后的列表
def rm_dup_none(l):
    return list(filter(None, sorted(set(l), key=l.index)))


# 定义一个名为 'filter_none' 的函数，该函数接受一个列表 'l' 作为参数
# 函数使用 'filter' 函数删除列表中所有的 'None' 元素
# 最后，函数返回处理后的列表
def filter_none(l):
    return list(filter(None, l))



# 定义一个名为GREPoreSeq的类
class GREPoreSeq:
    # 定义初始化函数，当类实例化时会被调用
    # 此处没有为类定义任何初始化属性
    def __init__(self):
        pass

    # 定义set_dir方法，用于创建并设置一个文件夹（目录）
    # 输入参数是文件夹名
    def set_dir(self, dir_name):
        """Setting directory"""
        # 检查这个文件夹是否已经存在
        if dir_name not in os.listdir():
            # 如果不存在，则创建该文件夹
            os.mkdir(dir_name)
            # 返回创建的文件夹名称
            return dir_name
        else:
            # 如果文件夹已经存在，则直接返回文件夹名称
            return dir_name

    # 定义seq_rev_com方法，用于获取输入序列的反向互补序列
    # 输入参数是序列信息，包括序列ID，序列本身，质量分数等
    def seq_rev_com(self, read):
        """reverse complement"""
        # 定义一个翻译表，用于将DNA序列转化为其反向互补序列????????
        trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
        # 获取反向互补序列
        rev_com_seq = read[1].translate(trantab)[::-1]
        # 生成新的序列ID，后缀加上'_R'表示这是反向互补序列
        rev_com_id = f'{read[0].strip().split()[0]}_R'
        # 生成新的序列信息
        linel = [rev_com_id, rev_com_seq, read[2], read[3][::-1]]
        # 返回新的序列信息
        return linel

    # 定义rev_com方法，用于获取输入序列的反向互补序列
    # 输入参数是一个序列
    def rev_com(self, seq):
        # 检查输入序列是否为空
        if seq:
            # 定义一个翻译表，用于将DNA序列转化为其反向互补序列
            trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
            # 返回反向互补序列
            return seq.translate(trantab)[::-1]
        else:
            # 如果输入序列为空，则返回None
            return

    # 定义input_info方法，用于从YAML文件中提取信息
    # 输入参数是两个YAML文件的路径
    def input_info(self, ref_yaml=None, demulti_yaml=None):
        """Extraction of excel information"""
        # 如果提供了ref_yaml文件的路径，则从该文件中提取信息
        if ref_yaml:
            self.ref_info = ri.r_input(ref_yaml)
        # 如果提供了demulti_yaml文件的路径，则从该文件中提取信息
        if demulti_yaml:
            self.demul_info = ri.r_input(demulti_yaml)
    # 定义mk_reference方法，用于生成参考序列文件
    def mk_reference(self):
        """Check and create the output folder"""
        # 记录开始创建FASTA文件的日志
        logger.info('Making FASTA file')
        # 创建名为'Reference'的文件夹，并保存其路径
        self.ref_dir = self.set_dir('Reference')

        # 创建一个空字典，用于存放参考序列ID和文件路径
        self.id_references = {}
        # 遍历所有的参考序列信息
        for ref_id, reference_info in self.ref_info.items():
            # 格式化参考序列ID
            reference_id = self.reform_id(ref_id)
            # 格式化参考序列
            seq = self.reform_seq(reference_info['sequence'])
            # 创建FASTA格式的字符串
            line = f'>{reference_id}\n'
            line += textwrap.fill(seq, 60)
            # 生成输出文件的路径
            fa_out = os.path.join(self.ref_dir, f'{reference_id}.fasta')
            # 将参考序列ID和输出文件路径保存到字典中
            self.id_references[reference_id] = fa_out
            # 打开输出文件，并将FASTA格式的字符串写入文件
            with open(fa_out, 'w') as f:
                f.write(line)

    # 定义mk_grepseq方法，用于生成grep序列，经过一定处理得到的用于搜索和比对的短序列
    def mk_grepseq(self, ref_seq, walk, step, seq_name='grepseq'):
        # 检查参考序列是否存在
        if ref_seq:
            # 创建名为'grepseqs'的文件夹，并保存其路径
            outpath = self.set_dir('grepseqs')
            # 获取参考序列的长度
            s, length = 0, len(ref_seq)
            # 根据步长和窗口大小，生成grep序列
            greps = [ref_seq[s + step * i:walk + step * i] for i in range(round(length / step))]
            # 从grep序列中筛选出长度等于窗口大小的序列
            grepws = [seq for seq in greps if len(seq) == walk]
            # 为筛选出的序列添加换行符
            grepseqws_lines = [seq+'\n' for seq in grepws]
            # 生成输出文件的路径
            grepseq_out = os.path.join(outpath, f"{seq_name}.txt")
            # 打开输出文件，并将筛选出的序列写入文件
            with open(grepseq_out, 'w') as g_out:
                g_out.writelines(grepseqws_lines)
            # 返回筛选出的序列
            return grepws
        else:
            # 如果参考序列不存在，则返回None
            return

    # 定义osearch方法，用于在目标序列中搜索子序列
    def osearch(self, seqs_used, seq_searched, match):
        # 初始化匹配数为0
        m = 0
        # 检查子序列列表是否存在
        if seqs_used:
            # 遍历子序列列表
            for seq in seqs_used:
                # 在目标序列中搜索子序列
                if re.search(seq, seq_searched):
                    # 如果找到子序列，则匹配数加1
                    m += 1
                    # 如果匹配数大于等于给定的阈值，则返回匹配数
                    if m >= match:
                        return m
        elif not seqs_used:
            # 如果子序列列表不存在，则返回True
            return True

    # 定义reform_seq方法，用于格式化序列
    def reform_seq(self, seq):
        # 检查序列是否存在
        if seq:
            # 去除序列中的空格、换行符和制表符，并将序列转化为大写
            return seq.replace(' ', '').replace('\n', '').replace('\t', '').upper()
        else:
            # 如果序列不存在，则返回None
            return

    # 定义reform_num方法，用于格式化数字
    def reform_num(self, num):
        # 检查数字是否存在
        if num:
            # 将数字转化为整数
            return int(num)
        else:
            # 如果数字不存在，则返回None
            return

    # 定义reform_id方法，用于格式化ID
    def reform_id(self, name_id):
        # 检查ID是否存在
        if name_id:
            # 将ID中的空格、句点和制表符替换为下划线
            return name_id.replace(' ', '_').replace('.', '_').replace('\t', '_')
        else:
            # 如果ID不存在，则返回None
            return

    # 定义disassemble方法，用于解析输入信息
    def disassemble(self):
        # 记录开始解析输入信息的日志
        logger.info('Disassembling input info')
        # 初始化各种属性
        self.fastq_ids = []
        self.id_ref, self.left_seqs, self.right_seqs, self.BCprimer_Fs, self.BClen_Fs, self.id_BCprimer_Rs, self.id_uniseq, self.id_BClen_R = {}, {}, {}, {}, {}, {}, {}, {}
        # 遍历所有的demulti信息
        for demuti_id, demuti_info in self.demul_info.items():
            # 格式化demulti ID
            demuti_id = self.reform_id(demuti_id)
            # 将demulti ID添加到列表中
            self.fastq_ids.append(demuti_id)
            # 解析并保存demulti信息
            self.left_seqs[f'{demuti_id}'] = self.reform_seq(demuti_info['left_seqs'])
            self.right_seqs[f'{demuti_id}'] = self.reform_seq(demuti_info['right_seqs'])
            self.BCprimer_Fs[f'{demuti_id}'] = self.reform_seq(demuti_info['BCprimer_F'])
            self.BClen_Fs[f'{demuti_id}'] = self.reform_num(demuti_info['BClen_F'])
            self.id_ref[f'{demuti_id}'] = self.reform_id(demuti_info['reference_id'])
            self.id_BCprimer_Rs[f'{demuti_id}'] = self.rev_com(self.reform_seq(demuti_info['BCprimer_R']))
            self.id_BClen_R[f'{demuti_id}'] = self.reform_num(demuti_info['BClen_R'])
            self.id_uniseq[f'{demuti_id}'] = self.reform_seq(demuti_info['unique_sequence'])
    # 准备grep序列
    def prepare_grepseq(self):
        # 记录日志：正在准备 grep 序列
        logger.info('Preparing Grepseqs')

        # 初始化两个字典，用于存储与各个 id 关联的 reads 和 grep 序列
        self.id_reads_dic, self.id_grepseq_dic = {}, {}

        # 遍历每一个 fastq 文件的 id
        for f_id in self.fastq_ids:
            # 从相应的字典中提取出各个 id 对应的序列和信息
            self.l_seq = self.left_seqs[f_id]
            self.r_seq = self.right_seqs[f_id]
            BCprimer_F = self.BCprimer_Fs[f_id]
            BCprimerR = self.id_BCprimer_Rs[f_id]
            BClen_F = self.BClen_Fs[f_id]
            BClen_R = self.id_BClen_R[f_id]
            uniseq = self.id_uniseq[f_id]

            # 如果 BClen_F 或 BClen_R 不存在，对应的 BCrange 设置为 -1，否则进行计算
            BCrange_F = -1 if not BClen_F else BClen_F - 4 + 11
            BCrange_R = -1 if not BClen_R else BClen_R - 4 + 11

            # 初始化每个 id 对应的 reads 列表
            self.id_reads_dic[f_id] = []

            # 根据是否存在对应的序列，创建相应的 grep 序列，并添加到 id_grepseq_dic 字典中
            # 如果某个序列不存在，对应的 grep 序列设置为 None
            # 生成 grep 序列的方法是 mk_grepseq，其中的参数有：
            # ref_seq: 参考序列
            # walk: 滑动窗口的长度
            # step: 滑动步长
            # seq_name: 用于生成 grep 序列文件的文件名
            if self.l_seq:
                self.id_grepseq_dic[f'{f_id}leftseq'] = self.mk_grepseq(ref_seq=self.l_seq, walk=20, step=7, seq_name=f'{f_id}leftseq')
            else:
                self.id_grepseq_dic[f'{f_id}leftseq'] = None
            if self.r_seq:
                self.id_grepseq_dic[f'{f_id}rightseq'] = self.mk_grepseq(ref_seq=self.r_seq, walk=20, step=7, seq_name=f'{f_id}rightseq')
            else:
                self.id_grepseq_dic[f'{f_id}rightseq'] = None
            if BCprimer_F:
                self.id_grepseq_dic[f'{f_id}BCprimer_F'] = self.mk_grepseq(ref_seq=BCprimer_F[:BCrange_F], walk=11,
                                                                           step=1, seq_name=f'{f_id}BCprimer_F')
            else:
                self.id_grepseq_dic[f'{f_id}BCprimer_F'] = None
            if BCprimerR:
                self.id_grepseq_dic[f'{f_id}BCprimer_R'] = self.mk_grepseq(ref_seq=BCprimerR[:BCrange_R], walk=11,
                                                                           step=1, seq_name=f'{f_id}BCprimer_R')
            else:
                self.id_grepseq_dic[f'{f_id}BCprimer_R'] = None
            if uniseq:
                steps = len(uniseq) // 20
                if steps <= 0:
                    steps = 1
                self.id_grepseq_dic[f'{f_id}uniseq'] = self.mk_grepseq(ref_seq=uniseq, walk=17, step=steps, seq_name=f'{f_id}uniseq')
            else:
                self.id_grepseq_dic[f'{f_id}uniseq'] = None

    def is_complete(self, f_id, read, read_rev_com, match):
        # 从 id_grepseq_dic 字典中提取各类 grep 序列
        gerpseqs_leftseq = self.id_grepseq_dic[f'{f_id}leftseq']
        gerpseqs_rightseq = self.id_grepseq_dic[f'{f_id}rightseq']
        grepseqs_BCprimerF = self.id_grepseq_dic[f'{f_id}BCprimer_F']
        grepseqs_BCprimerR = self.id_grepseq_dic[f'{f_id}BCprimer_R']
        grepseqs_uniseq = self.id_grepseq_dic[f'{f_id}uniseq']

        # 初始化左右端的范围
        left_range = 50
        right_range = 50

        # 如果左端序列存在，则其范围设置为该序列的长度
        if self.l_seq:
            left_range = len(self.l_seq)

        # 如果右端序列存在，则其范围设置为该序列的长度
        if self.r_seq:
            right_range = len(self.r_seq)

        # 通过 osearch 方法在左端 grep 序列中搜索目标序列或其反向互补序列
        match_nf_leftseq = self.osearch(seqs_used=gerpseqs_leftseq, seq_searched=read[1][:left_range], match=match)
        match_nr_leftseq = self.osearch(seqs_used=gerpseqs_leftseq, seq_searched=read_rev_com[1][:left_range],
                                        match=match)

        # 如果左端 grep 序列不存在，而右端 grep 序列存在
        if not gerpseqs_leftseq and gerpseqs_rightseq:
            match_nr_rightseq_ = self.osearch(seqs_used=gerpseqs_rightseq, seq_searched=read_rev_com[1][-right_range:],
                                              match=match)

            # 如果在右端 grep 序列中找到了匹配项，进一步搜索 BCprimerF、BCprimerR 和 uniseq 序列
            if match_nr_rightseq_:
                match_n_BCprimerF = self.osearch(seqs_used=grepseqs_BCprimerF, seq_searched=read_rev_com[1][:21],
                                                 match=match)
                if match_n_BCprimerF:
                    match_n_BCprimerR = self.osearch(seqs_used=grepseqs_BCprimerR, seq_searched=read_rev_com[1][-20:],
                                                     match=match)
                    if match_n_BCprimerR:
                        match_n_uniseq = self.osearch(seqs_used=grepseqs_uniseq,
                                                      seq_searched=read_rev_com[1][left_range:-right_range],
                                                      match=match)
                        # 如果在所有的序列中都找到了匹配项，返回"NR"
                        if match_n_uniseq:
                            return "NR"

        # 如果在左端 grep 序列中找到了匹配项，进一步搜索右端 grep 序列和 BCprimerF、BCprimerR、uniseq 序列
        if match_nf_leftseq:
            match_nf_rightseq = self.osearch(seqs_used=gerpseqs_rightseq, seq_searched=read[1][-right_range:],
                                             match=match)
            if match_nf_rightseq:
                match_nf_BCprimerF = self.osearch(seqs_used=grepseqs_BCprimerF, seq_searched=read[1][:21], match=match)
                if match_nf_BCprimerF:
                    match_nf_BCprimerR = self.osearch(seqs_used=grepseqs_BCprimerR, seq_searched=read[1][-20:],
                                                      match=match)
                    if match_nf_BCprimerR:
                        match_nf_uniseq = self.osearch(seqs_used=grepseqs_uniseq,
                                                       seq_searched=read[1][left_range:-right_range], match=match)
                        # 如果在所有的序列中都找到了匹配项，返回"NF"
                        if match_nf_uniseq:
                            return "NF"
        # 如果在右端 grep 序列中找到了匹配项，
        elif match_nr_leftseq:
            match_nr_rightseq = self.osearch(seqs_used=gerpseqs_rightseq, seq_searched=read_rev_com[1][-right_range:],
                                             match=match)
            if match_nr_rightseq:
                match_nr_BCprimerF = self.osearch(seqs_used=grepseqs_BCprimerF, seq_searched=read_rev_com[1][:21],
                                                  match=match)
                if match_nr_BCprimerF:
                    match_nr_BCprimerR = self.osearch(seqs_used=grepseqs_BCprimerR,
                                                      seq_searched=read_rev_com[1][-20:],
                                                      match=match)
                    if match_nr_BCprimerR:
                        match_nr_uniseq = self.osearch(seqs_used=grepseqs_uniseq,
                                                       seq_searched=read_rev_com[1][left_range:-right_range],
                                                       match=match)
                        if match_nr_uniseq:
                            return "NR"
#以上代码主要实现了 demultiplexing，写入 fastq 文件以及进行可视化操作的功能。demultiplexing 通过检查序列是否完整并分类。写入 fastq 文件则是将 demultiplex 后的结果写入到相应的 fastq 文件中，并对结果进行简单的统计。最后，可视化函数则是通过一个叫 vis 的模块进行可视化操作。
    # 分流
    def demultiplex(self, fastq_file, match=1):
        # 打印开始进行demultiplex的日志信息
        logger.info(f'Demultiplexing {fastq_file}')

        # 初始化计数器
        nt = 0
        # 记录开始时间
        start = time.time()
        # 对 fastq 文件中的每一条序列进行操作
        for read in read_fastq(fastq_file):
            # 获取反向互补序列
            revcom_line = self.seq_rev_com(read)

            # 对每一个 fastq id 进行操作
            for f_id in self.fastq_ids:
                # 检查该序列是否完整
                result = self.is_complete(f_id, read, revcom_line, match)

                # 如果结果为 "NF"，将该序列添加到对应 id 的列表中
                if result == "NF":
                    line = f"{read[0]}\n{read[1]}\n{read[2]}\n{read[3]}\n"
                    self.id_reads_dic[f_id].append(line)

                # 如果结果为 "NR"，将该反向互补序列添加到对应 id 的列表中
                elif result == "NR":
                    line_rc = f"{revcom_line[0]}\n{revcom_line[1]}\n{revcom_line[2]}\n{revcom_line[3]}\n"
                    self.id_reads_dic[f_id].append(line_rc)

            # 更新计数器
            nt += 1

            # 每处理100000条序列，打印一次进度信息
            if nt % 100000 == 0:
                logger.info("Processed %d reads in %.1f minutes" % (nt, (time.time() - start) / 60))

    def write_fastq(self):
        # 设置输出文件夹
        self.demulti_outpath = self.set_dir('Demultiplexed')
        self.demultis, stats = [], []
        self.demulti_random200s = {}
        stats = ['file\treads\n']
        # 对 id_reads_dic 中每个 id 对应的序列进行操作
        for id, reads in self.id_reads_dic.items():
            # 定义输出文件的路径
            out = os.path.join(self.demulti_outpath, id + '.fastq')
            out_random200 = os.path.join(self.demulti_outpath, id + '_random200.fastq')
            self.demultis.append(out)
            self.demulti_random200s[out] = None
            # 打开文件并写入序列
            with open(out, 'w') as f_out:
                if reads:
                    lines = reads[:-2]
                    lines.append(reads[-1].strip())
                    f_out.writelines(lines)
                    # 如果序列数量大于200，随机取200条写入文件
                    if len(lines) > 200:
                        cmd_200 = f"seqkit sample -n 200 {out} -o {out_random200}"
                        self.demulti_random200s[out] = out_random200
                        os.system(cmd_200)
                else:
                    continue
            # 打印信息，告知序列已经写入文件
            logger.info(f"{len(reads)} reads are written to the {out}")
            stats.append(f"{out}\t{len(reads)}\n")
            # 运行seqkit命令，生成统计图
            seqkit_watch_cmd = f'seqkit watch -Q --fields ReadLen {out} -O {out.strip(".fastq")}.pdf'
            subprocess.call(seqkit_watch_cmd, shell=True)
        # 将统计信息写入文件
        stats_out = os.path.join(self.demulti_outpath, 'Demultiplex_stats.txt')
        with open(stats_out, 'w') as f_stat:
            f_stat.writelines(stats)

    def write_to_file(self, file, lines):
        for line in lines:
            file.write(line)

    def demul_P_N(self):
        print("=======================demul_P_N开始运行============================")

               # 获取当前脚本所在的目录
        script_dir = os.path.dirname(os.path.abspath(__file__))
        print(f"========={script_dir}=========")
       # 获取上级目录
        parent_dir = os.path.dirname(script_dir)
        # 构建 Demultiplexed 文件夹的路径
        demultiplexed_dir = os.path.join(parent_dir, "Demultiplexed")
        # 检查 Demultiplexed 文件夹是否存在
        fastq_files = []
        if os.path.exists(demultiplexed_dir):
            # 遍历 Demultiplexed 文件夹中的所有文件
            for file_name in os.listdir(demultiplexed_dir):
                # 检查文件是否以 .fastq 结尾
                if file_name.endswith(".fastq"):
                    # 构建 fastq 文件的完整路径
                    fastq_files.append(os.path.join(demultiplexed_dir, file_name))
        else:
            print("==========fastq_files为空================")

        for file_name in fastq_files:
            # 提取文件名的基本部分，以便在输出文件名中使用
            base_file_name = os.path.splitext(os.path.basename(file_name))[0]

            with open(file_name, 'r') as original_file, \
                    open(os.path.join(demultiplexed_dir, f'{base_file_name}-N.fastq'), 'w') as R_file, \
                    open(os.path.join(demultiplexed_dir, f'{base_file_name}-P.fastq'), 'w') as no_R_file:

                lines_for_sequence = []
                for line in original_file:
                    lines_for_sequence.append(line)

                    # 检查我们是否看到了这个序列的所有4行
                    if len(lines_for_sequence) == 4:
                        # 如果这个序列以_R结尾
                        if lines_for_sequence[0].strip().endswith('_R'):
                            self.write_to_file(R_file, lines_for_sequence)
                        else:
                            self.write_to_file(no_R_file, lines_for_sequence)

                        # 为下一个序列重置
                        lines_for_sequence = []

      
    def visualization(self):
        # 设置输出目录
        self.visua_abspath = self.set_dir('Visualization')
        sorted_bams = []
        # 对每一个demultiplexed序列进行可视化操作
        for demulti_id, demulti in zip(self.fastq_ids, self.demultis):
            ref_id = self.id_ref[demulti_id]
            reference = self.id_references[ref_id]
            demulti_random200 = self.demulti_random200s[demulti]
            # 进行可视化操作，返回结果
            sorted_bam = vis.visualizing(reference, demulti, self.visua_abspath, demulti_random200)
            sorted_bams.append(sorted_bam)
        # 打印完成信息
        logger.info("Visualization done")


def parse_args():
    # 创建argparse解析器
    parser = argparse.ArgumentParser()
    # 添加子解析器
    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')
    # 添加‘all’子解析器
    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('-d', '--demultiplexinfo', help="Specify the DemultiplexInfo file", required=True)
    all_parser.add_argument('-n', '--nanopore', help="Specify the Nanopore sequencing data file", required=True)
    all_parser.add_argument('-r', '--referenceinfo', help="Specify the ReferenceInfo file", required=True)

    # 添加‘mkreference’子解析器
    mkreference_parser = subparsers.add_parser('mkreference',
                                              help='Create FASTA files based on the information in the "ReferenceInfo"')
    mkreference_parser.add_argument('-r', '--referenceinfo', help="Specify the ReferenceInfo file", required=True)

    # 添加‘mkgrepseq’子解析器
    mkgrepseq_parser = subparsers.add_parser('mkgrepseq', help='Create GREPseq files')
    mkgrepseq_parser.add_argument('-d', '--demultiplexinfo', help="Specify the DemultiplexInfo file", required=True)

    # 添加‘demultiplex’子解析器
    demultiplex_parser = subparsers.add_parser('demultiplex',
                                               help='Demultiplex FASTQ files based on the information in the "DemultiplexInfo"')
    demultiplex_parser.add_argument('-n', '--nanopore', help="Specify the Nanopore sequencing data file", required=True)
    demultiplex_parser.add_argument('-d', '--demultiplexinfo', help="Specify the DemultiplexInfo file", required=True)

    # 添加‘visualize’子解析器
    visualize_parser = subparsers.add_parser('visualize', help='Minimap2 align FASTQ files, '
                                                               'Samtools generate the sorted.bam and bai files '
                                                               'that IGV needed')
    visualize_parser.add_argument('-a', '--fasta', help='Specify the reference FASTA file', required=True)
    visualize_parser.add_argument('-q', '--fastq', help='Specify demultiplexed FASTQ files', required=True, nargs='+')

    return parser.parse_args()  # 解析参数并返回

def main():
    def if_list(var):
        # 如果变量是列表，直接返回；否则，将其放入一个列表中并返回
        if isinstance(var, list):
            return var
        else:
            return [var]

    args = parse_args()  # 解析命令行参数
    if args.command == 'all':
        # 运行pipeline的所有步骤
        g = GREPoreSeq()
        # 从YAML文件中提取信息
        g.input_info(ref_yaml=args.referenceinfo, demulti_yaml=args.demultiplexinfo)
        g.mk_reference() #生成参考序列文件
        g.disassemble() #用于解析输入信息
        g.prepare_grepseq()  #准备grep序列
        g.demultiplex(fastq_file=args.nanopore, match=1) #分流
        g.write_fastq() #写fastq
        g.demul_P_N()
        g.visualization() #可视化

    elif args.command == 'mkreference':
        """
        Run just mkreference step given the RefInfo
        """
        # 只运行mkreference步骤
        g = GREPoreSeq()
        g.input_info(ref_yaml=args.referenceinfo)
        g.mk_reference()

    elif args.command == 'mkgrepseq':
        # 运行mkgrepseq步骤
        g = GREPoreSeq()
        g.input_info(demulti_yaml=args.demultiplexinfo)
        g.disassemble()
        g.prepare_grepseq()

    elif args.command == 'demultiplex':
        # 运行demultiplex步骤
        g = GREPoreSeq()
        g.input_info(demulti_yaml=args.demultiplexinfo)
        g.disassemble()
        g.prepare_grepseq()
        g.demultiplex(fastq_file=args.nanopore, match=1)
        g.write_fastq()
        g.demul_P_N()

    elif args.command == 'visualize':
        # 运行visualize步骤
        g = GREPoreSeq()
        outpath = g.set_dir('Visualization')
        fastqs = if_list(args.fastq)
        for demuti in fastqs:
            vis.visualizing(args.fasta, demuti, outpath)

#以上的代码创建了一个命令行界面，通过命令行参数允许用户指定需要运行的步骤以及相关参数。
# 其中包含了all、mkreference、mkgrepseq、demultiplex和visualize等五个命令，
# 分别对应着管道中的不同步骤。在main函数中，这些命令会被相应地执行。
if __name__ == '__main__':
    main()
