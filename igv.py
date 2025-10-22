# import os
# import sys

# # 设置默认编码为 UTF-8
# os.environ["PYTHONIOENCODING"] = "utf-8"
# sys.stdout.reconfigure(encoding='utf-8')

# # def windows_to_wsl_path(windows_path):
# #     """Convert a Windows path to a WSL path."""
# #     print(f'windows_path:{windows_path}')
# #     windows_path = windows_path.replace("\\", "/")
# #     drive, rest_of_path = windows_path.split(":/", 1)
# #     wsl_path = f"/mnt/{drive.lower()}/{rest_of_path}"
# #     return wsl_path

# # Define folder paths
# current_path = os.getcwd()
# print(current_path)
# # 获取 reference 文件夹的路径

# reference_folder = os.path.join(current_path, "Reference")  # Replace with the folder containing reference files
# visualization_folder = os.path.join(current_path, "Visualization") # Replace with the folder containing BAM files
# reference_folder2 = os.path.join(current_path, "Demultiplexed")
# visualization_folder2 = os.path.join(current_path, "Demultiplexed/Visualization")
# output_dir = os.path.join(current_path, "IGV")  # Replace with the output folder path
# batch_file = os.path.join(current_path, "igv.txt") # Replace with your batch file path

# # Automatically detect files
# fasta_files = [f for f in os.listdir(reference_folder) if f.endswith(".fasta")]
# fasta_files2 = [f for f in os.listdir(reference_folder2) if f.endswith("RC.fasta")]
# bam_files = [f for f in os.listdir(visualization_folder) if f.endswith(".sorted.bam")]
# bam_files2 = [f for f in os.listdir(visualization_folder2) if f.endswith("-L_200-RC.sorted.bam")]

# # Ensure there is at least one reference genome and BAM file
# if len(fasta_files) == 0:
#     raise ValueError("No .fasta files found in the reference folder!")
# if len(bam_files) == 0:
#     raise ValueError("No .sorted.bam files found in the visualization folder!")
# if len(fasta_files2) == 0:
#     raise ValueError("No .fasta files found in the Demultiplexed folder!")
# if len(bam_files2) == 0:
#     raise ValueError("No .sorted.bam files found in the Demultiplexed/visualization folder!")


# # Convert paths to WSL format
# # reference_folder_wsl = windows_to_wsl_path(reference_folder)
# # visualization_folder_wsl = windows_to_wsl_path(visualization_folder)
# # reference_folder_wsl2 = windows_to_wsl_path(reference_folder2)
# # visualization_folder_wsl2 = windows_to_wsl_path(visualization_folder2)
# # output_dir_wsl = windows_to_wsl_path(output_dir)
# # batch_file_wsl = windows_to_wsl_path(batch_file)

# reference_folder_wsl = reference_folder
# visualization_folder_wsl = visualization_folder
# reference_folder_wsl2 = reference_folder2
# visualization_folder_wsl2 = visualization_folder2
# output_dir_wsl = output_dir
# batch_file_wsl = batch_file

# # Create a mapping of reference files to BAM files
# # Assuming the BAM files are named with the same prefix as the reference files
# reference_to_bam = {}
# for fasta_file in fasta_files:
#     prefix = os.path.splitext(fasta_file)[0]  # Get the file name without extension
#     matched_bam_files = [
#         bam for bam in bam_files if bam.startswith(prefix) and not bam.startswith(f"{prefix}-L")
#     ]
#     if matched_bam_files:
#         reference_to_bam[fasta_file] = matched_bam_files
#     else:
#         print(f"Warning: No BAM files found for reference {fasta_file}")
# reference_to_bam2 = {}
# for fasta_file in fasta_files2:
#     prefix = os.path.splitext(fasta_file)[0][:5]  # Get the file name without extension
#     # print(prefix)
#     matched_bam_files = [bam for bam in bam_files2 if bam.startswith(prefix)]
#     if matched_bam_files:
#         reference_to_bam2[fasta_file] = matched_bam_files
#     else:
#         print(f"Warning: No BAM files found for reference {fasta_file}")


# # Create batch file content
# with open(batch_file, "w", encoding="utf-8") as f:
#     f.write("new\n")
#     f.write(f"snapshotDirectory {output_dir_wsl}\n\n")

#     for reference, bam_list in reference_to_bam.items():
#         reference_path = os.path.join(reference_folder, reference)
#         reference_path_wsl = reference_path

#         f.write(f"genome {reference_path_wsl}\n")
#         for bam_file in bam_list:
#             bam_path = os.path.join(visualization_folder, bam_file)
#             bam_path_wsl = bam_path
#             bam_name = os.path.basename(bam_file).replace(".bam", "")

#             f.write("clear\n")
#             f.write(f"load {bam_path_wsl}\n")
#             f.write("colorBy READ_STRAND\n")
#             f.write("squish\n")  # Set tracks to squished display
#             f.write("maxPanelHeight 800\n")  # Set image height
#             f.write(f"snapshot {bam_name}_by_strand.png\n\n")

#     for reference, bam_list in reference_to_bam2.items():
#         reference_path = os.path.join(reference_folder2, reference)
#         reference_path_wsl = reference_path

#         f.write(f"genome {reference_path_wsl}\n")
#         for bam_file in bam_list:
#             bam_path = os.path.join(visualization_folder2, bam_file)
#             bam_path_wsl = bam_path
#             bam_name = os.path.basename(bam_file).replace(".bam", "")

#             f.write("clear\n")
#             f.write(f"load {bam_path_wsl}\n")
#             f.write("colorBy READ_STRAND\n")
#             f.write("squish\n")  # Set tracks to squished display
#             f.write("maxPanelHeight 800\n")  # Set image height
#             f.write(f"snapshot {bam_name}.png\n\n")

#     f.write("exit\n")

# print(f"Batch file generated: {batch_file}")
# print(f"Batch file (WSL path): {batch_file_wsl}")

import os
import sys

# 设置默认编码为 UTF-8
os.environ["PYTHONIOENCODING"] = "utf-8"
sys.stdout.reconfigure(encoding="utf-8")

# 当前路径
current_path = os.getcwd()
print(f"当前工作目录: {current_path}")

# 文件夹路径
reference_folder     = os.path.join(current_path, "Reference")
visualization_folder = os.path.join(current_path, "Visualization")

reference_folder2     = os.path.join(current_path, "Demultiplexed")
visualization_folder2 = os.path.join(current_path, "Demultiplexed/Visualization")

output_dir  = os.path.join(current_path, "IGV")
batch_file  = os.path.join(current_path, "igv.txt")

# 自动检测 Reference/Visualization
fasta_files = [f for f in os.listdir(reference_folder) if f.endswith(".fasta")]
bam_files   = [f for f in os.listdir(visualization_folder) if f.endswith(".sorted.bam")]

if len(fasta_files) == 0:
    raise ValueError("❌ 未在 Reference 文件夹中找到 .fasta 文件！")
if len(bam_files) == 0:
    raise ValueError("❌ 未在 Visualization 文件夹中找到 .sorted.bam 文件！")

# 自动检测 Demultiplexed（可选）
fasta_files2, bam_files2 = [], []
if os.path.isdir(reference_folder2) and os.path.isdir(visualization_folder2):
    fasta_files2 = [f for f in os.listdir(reference_folder2) if f.endswith("RC.fasta")]
    bam_files2   = [f for f in os.listdir(visualization_folder2) if f.endswith("-L_200-RC.sorted.bam")]
    if len(fasta_files2) == 0 or len(bam_files2) == 0:
        print("⚠️ Demultiplexed 文件夹存在，但未找到 RC.fasta 或 -L_200-RC.sorted.bam 文件，将跳过。")
else:
    print("⚠️ 未找到 Demultiplexed 相关目录，跳过这一部分。")

# 建立 fasta ↔ BAM 映射
reference_to_bam = {}
for fasta_file in fasta_files:
    prefix = os.path.splitext(fasta_file)[0]
    matched_bams = [bam for bam in bam_files if bam.startswith(prefix) and not bam.startswith(f"{prefix}-L")]
    if matched_bams:
        reference_to_bam[fasta_file] = matched_bams
    else:
        print(f"⚠️ Warning: 没找到 {fasta_file} 对应的 BAM 文件")

reference_to_bam2 = {}
for fasta_file in fasta_files2:
    prefix = os.path.splitext(fasta_file)[0][:5]  # 前缀匹配
    matched_bams = [bam for bam in bam_files2 if bam.startswith(prefix)]
    if matched_bams:
        reference_to_bam2[fasta_file] = matched_bams
    else:
        print(f"⚠️ Warning: 没找到 {fasta_file} 对应的 BAM 文件")

# 写入 IGV 批处理文件
os.makedirs(output_dir, exist_ok=True)
with open(batch_file, "w", encoding="utf-8") as f:
    f.write("new\n")
    f.write(f"snapshotDirectory {output_dir}\n\n")

    # 第一部分
    for reference, bam_list in reference_to_bam.items():
        reference_path = os.path.join(reference_folder, reference)
        f.write(f"genome {reference_path}\n")
        for bam_file in bam_list:
            bam_path = os.path.join(visualization_folder, bam_file)
            bam_name = os.path.basename(bam_file).replace(".bam", "")

            f.write("clear\n")
            f.write(f"load {bam_path}\n")
            f.write("colorBy READ_STRAND\n")
            f.write("squish\n")
            f.write("maxPanelHeight 800\n")
            f.write(f"snapshot {bam_name}_by_strand.png\n\n")

    # 第二部分（可选）
    if reference_to_bam2:
        for reference, bam_list in reference_to_bam2.items():
            reference_path = os.path.join(reference_folder2, reference)
            f.write(f"genome {reference_path}\n")
            for bam_file in bam_list:
                bam_path = os.path.join(visualization_folder2, bam_file)
                bam_name = os.path.basename(bam_file).replace(".bam", "")

                f.write("clear\n")
                f.write(f"load {bam_path}\n")
                f.write("colorBy READ_STRAND\n")
                f.write("squish\n")
                f.write("maxPanelHeight 800\n")
                f.write(f"snapshot {bam_name}.png\n\n")

    f.write("exit\n")

print(f"✅ IGV 批处理文件生成完成: {batch_file}")
print(f"📂 截图将保存在: {output_dir}")
