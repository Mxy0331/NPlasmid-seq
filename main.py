# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

import subprocess
import os
import shutil
import sys
import glob
import matplotlib.pyplot as plt
import re
import csv
import matplotlib
matplotlib.use("Agg")   # 必须放在 import matplotlib.pyplot 之前
import matplotlib.pyplot as plt

# def run_script(script_path):
#     """运行一个脚本"""
#     print(f"正在运行: {script_path}")
#     try:
#         subprocess.run([sys.executable, script_path], check=True)
#     except subprocess.CalledProcessError as e:
#         print(f"运行 {script_path} 出错: {e}")
#         sys.exit(1)
def run_script(script_path):
    """运行一个脚本"""
    print(f"正在运行: {script_path}")
    try:
        script_dir = os.path.dirname(script_path)
        subprocess.run([sys.executable, script_path], check=True, cwd=script_dir)
    except subprocess.CalledProcessError as e:
        print(f"运行 {script_path} 出错: {e}")
        sys.exit(1)


def run_command(cmd_list, use_shell=False):
    """运行命令行"""
    print(f"正在运行命令: {' '.join(cmd_list)}")
    try:
        subprocess.run(cmd_list, check=True, shell=use_shell)
    except subprocess.CalledProcessError as e:
        print(f"命令执行出错: {e}")
        sys.exit(1)

def copy_if_exists(src, dst):
    """复制文件，如果存在"""
    if os.path.exists(src):
        shutil.copy(src, dst)
        print(f"已复制: {src} -> {dst}")
    else:
        print(f"警告: {src} 未生成")


def count_fastq_reads(fq_file):
    """统计 fastq 文件中的 reads 数"""
    with open(fq_file, "r") as f:
        for i, _ in enumerate(f):
            pass
        return (i + 1) // 4 if 'i' in locals() else 0

import matplotlib.pyplot as plt
import re

def count_fastq_reads(fq_file):
    """统计 fastq 文件中的 reads 数"""
    with open(fq_file, "r") as f:
        for i, _ in enumerate(f):
            pass
        return (i + 1) // 4 if 'i' in locals() else 0

def plot_contam_clusters(consen_dir, data_dir, result_dir, threshold=0.005,
                         title_fontsize=25, label_fontsize=25, pct_fontsize=25):
    os.makedirs(result_dir, exist_ok=True)  # 确保目录存在
    """绘制每个样本的污染聚类饼图"""
    fastq_files = glob.glob(os.path.join(consen_dir, "*.fastq"))
    if not fastq_files:
        print("[WARNING] No fastq found in consen/, skip plotting.")
        return

    sample_pattern = re.compile(r"(P\d{4})")
    sample_groups = {}
    for fq in fastq_files:
        fname = os.path.basename(fq)
        m = sample_pattern.search(fname)
        if not m:
            print(f"[SKIP] {fname} has no sample ID")
            continue
        sample_id = m.group(1)
        sample_groups.setdefault(sample_id, []).append(fq)

    for sample_id, files in sample_groups.items():
        # 找原始 data fastq
        raw_fq = os.path.join(data_dir, f"{sample_id}.fastq")
        if not os.path.exists(raw_fq):
            print(f"[WARNING] Missing {raw_fq}, skip {sample_id}")
            continue

        total_reads = count_fastq_reads(raw_fq)
        if total_reads == 0:
            print(f"[WARNING] {sample_id} total reads = 0, skip")
            continue

        # 记录污染聚类
        contam_labels = []
        contam_sizes = []

        for idx, fq in enumerate(files, start=1):
            cluster_reads = count_fastq_reads(fq)
            proportion = cluster_reads / total_reads
            if proportion >= threshold:
                contam_labels.append(f"Cluster{idx} ({proportion*100:.2f}%)")
                contam_sizes.append(cluster_reads)
                print(f"[Contam] {sample_id} Cluster{idx}: {proportion*100:.2f}%")
            else:
                print(f"[Non-contam] {sample_id} Cluster{idx}: {proportion*100:.2f}% (ignored)")

        # 绘图
        plt.figure(figsize=(7, 7))
        if contam_sizes:  # 至少一个污染类
            colors = plt.cm.Set3(range(len(contam_sizes)))  # 柔和调色板
            wedges, texts, autotexts = plt.pie(
                contam_sizes,
                labels=contam_labels,
                colors=colors,
                autopct="%.2f%%",
                startangle=90,
                textprops={'fontsize': pct_fontsize}
            )
            plt.title(f"{sample_id} contamination clusters", fontsize=title_fontsize, fontweight="bold")
            for t in texts:
                t.set_fontsize(label_fontsize)
        else:  # 无污染
            plt.pie([1], labels=["No contamination"], colors=["#66c2a5"])
            plt.text(0, 0, f"{sample_id}\nNo contamination",
                     ha="center", va="center",
                     fontsize=title_fontsize, fontweight="bold")
            plt.title(f"{sample_id} contamination clusters", fontsize=title_fontsize, fontweight="bold")

        out_path = os.path.join(result_dir, f"{sample_id}_contam_piechart.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[OK] Saved {out_path}")

search_seq = ""   # 比如 "GCACTGCAGAGATGGATAACCA"，为空则走普通污染鉴别

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    word_dir = os.path.join(base_dir, "word")

    if not os.path.isdir(word_dir):
        print("错误: 找不到 word 文件夹")
        sys.exit(1)

    # # Step 1: 运行 CreateExcel.py
    # create_excel = os.path.join(word_dir, "CreateExcel.py")
    # run_script(create_excel)

    # # Step 2: 运行 genenrateYaml.py
    # generate_yaml = os.path.join(word_dir, "genenrateYaml.py")
    # run_script(generate_yaml)

    # # Step 3: 复制 yaml 文件
    # for yaml_name in ["refinfo.yaml", "Demuinfo.yaml"]:
    #     src = os.path.join(word_dir, yaml_name)
    #     dst = os.path.join(base_dir, yaml_name)
    #     copy_if_exists(src, dst)

    # # Step 4: 自动识别 fastq.gz 并运行 greporeseq
    # fastq_files = glob.glob(os.path.join(base_dir, "*.fastq.gz"))
    # if len(fastq_files) == 0:
    #     print("错误: 未找到任何 .fastq.gz 文件")
    #     sys.exit(1)
    # elif len(fastq_files) > 1:
    #     print("错误: 找到多个 .fastq.gz 文件，请确保只有一个")
    #     sys.exit(1)

    # fastq_file = os.path.basename(fastq_files[0])
    # greporeseq_script = os.path.join(base_dir, "greporeseq1_6", "greporeseq.py")

    # run_command([
    #     sys.executable, greporeseq_script,
    #     "all",
    #     "-n", fastq_file,
    #     "-d", "Demuinfo.yaml",
    #     "-r", "refinfo.yaml"
    # ])

    # # Step 5: 运行 igv.py 生成 igv.txt
    # igv_script = os.path.join(base_dir, "igv.py")
    # run_script(igv_script)

    # # Step 6: 调用 IGV 批处理执行截图
    # igv_batch_file = os.path.join(base_dir, "igv.txt")
    # if os.path.exists(igv_batch_file):
    #     run_command(["igv.sh", "-b", igv_batch_file])
    # else:
    #     print("警告: 未找到 igv.txt，跳过 IGV 执行")

    # # Step 7: 生成共识序列 (longread_umi racon)
    # consensus_dir = os.path.join(base_dir, "consensus")
    # consen_subdir = os.path.join(consensus_dir, "consen")
    # racon_outdir = os.path.join(consensus_dir, "racon")

    # os.makedirs(consen_subdir, exist_ok=True)
    # os.makedirs(racon_outdir, exist_ok=True)

    # # 将 Demultiplexed 下的 random200.fastq 复制到 consensus/consen
    # demux_dir = os.path.join(base_dir, "Demultiplexed")
    # fastq_sources = glob.glob(os.path.join(demux_dir, "*random200.fastq"))

    # i = 1
    # for src in fastq_sources:
    #     orig_name = os.path.basename(src)
    #     # 去掉后缀
    #     file_name, ext = os.path.splitext(orig_name)
    #     # 新名字格式：umi{i}_原文件名_bins.fastq
    #     new_name = f"umi{i}_{file_name}_bins{ext}"
    #     dst = os.path.join(consen_subdir, new_name)
    #     shutil.copy(src, dst)
    #     print(f"已复制并改名: {src} -> {dst}")
    #     i += 1

    # # 调用 longread_umi 生成 consensus_racon.fa
    # run_command([
    # "conda", "run", "-n", "longread_umi",
    # "longread_umi", "consensus_racon",
    # "-d", consen_subdir,
    # "-o", racon_outdir,
    # "-r", "3",
    # "-t", "16",
    # "-p", "map-ont"
    # ])

    # Step 8: 运行 BLAST 报告脚本
    blast_script = os.path.join(base_dir, "blast_result1_5_only_blastresults.py")
    if os.path.exists(blast_script):
        run_command([sys.executable, blast_script])
    else:
        print("警告: 未找到 blast_result1_5_only_blastresults.py，跳过 BLAST 报告")

   # Step 9: 污染鉴别
    contam_dir = os.path.join(base_dir, "contam_detect")

    if search_seq.strip():  
    # 如果提供了 search_seq，运行 20nt 模式
        contam20_dir = os.path.join(contam_dir, "20nt")
        contam20_result_dir = os.path.join(contam20_dir, "result")
        os.makedirs(contam20_result_dir, exist_ok=True)

        gridsearch_script = os.path.join(contam20_dir, "gridsearch_3params.py")
        if os.path.exists(gridsearch_script):
            # 修改 gridsearch_3params.py 中的 search_seq
            with open(gridsearch_script, "r") as f:
                lines = f.readlines()
            with open(gridsearch_script, "w") as f:
                for line in lines:
                    if line.strip().startswith("search_seq = "):
                        f.write(f'search_seq = "{search_seq}"\n')
                    else:
                        f.write(line)
            print(f"[INFO] 已将 gridsearch_3params.py 中的 search_seq 修改为 {search_seq}")

            # 运行脚本 —— 工作目录切换到 20nt
            cur_dir = os.getcwd()
            os.chdir(contam20_dir)
            try:
                run_command([sys.executable, gridsearch_script])
            finally:
                os.chdir(cur_dir)
        else:
            print("错误: 未找到 gridsearch_3params.py")

    else:
        # 如果没有 search_seq，走普通污染鉴别
        contam_data_dir = os.path.join(contam_dir, "data")
        contam_result_dir = os.path.join(contam_dir, "result")
        os.makedirs(contam_result_dir, exist_ok=True)

        # (a) 移动 Demultiplexed 中的原始 fastq (只保留 Pxxxx.fastq)
        demux_dir = os.path.join(base_dir, "Demultiplexed")
        pattern = re.compile(r"^P\d{4}\.fastq$")
        demux_fastqs = glob.glob(os.path.join(demux_dir, "*.fastq"))
        for fq in demux_fastqs:
            fname = os.path.basename(fq)
            if pattern.match(fname):
                dst = os.path.join(contam_data_dir, fname)
                shutil.copy(fq, dst)
                print(f"[OK] 已复制 {fq} -> {dst}")
            else:
                print(f"[SKIP] 跳过 {fname}")

        # (b) 移动 Reference 下的所有 .fasta
        ref_fastas = glob.glob(os.path.join(base_dir, "Reference", "*.fasta"))
        for fa in ref_fastas:
            dst = os.path.join(contam_data_dir, os.path.basename(fa))
            shutil.copy(fa, dst)
            print(f"[OK] 已复制 {fa} -> {dst}")

        # (c) 运行普通污染鉴别脚本
        contam_script = os.path.join(contam_dir, "scripts", "NPlasmid-seq.py")
        if os.path.exists(contam_script):
            cur_dir = os.getcwd()
            os.chdir(contam_result_dir)
            try:
                run_command([sys.executable, contam_script, "all", "-d", "../data", "-n", "5"])
            finally:
                os.chdir(cur_dir)
        else:
            print("警告: 未找到 NPlasmid-seq.py，跳过污染鉴别")

        # (d) 绘制饼图
        consen_dir = os.path.join(contam_result_dir, "consen")
        plot_contam_clusters(consen_dir, contam_data_dir, contam_result_dir, threshold=0.005)

    print("全部流程完成 ✅")

if __name__ == "__main__":
    main()
