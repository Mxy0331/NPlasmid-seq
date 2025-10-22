import os

print("将指定的文件类型的文件改名为统一格式！")
sub = input(f"输入文件后缀：")
files = list(os.listdir())
files.remove('Reformat.py')
i = 1
for file in files:
    if file.endswith(sub):
        file_name = file.strip(sub)
        file_rename = f'umi{i}_{file_name}_bins{sub}'
        os.rename(file, file_rename)
        i += 1

print('done')
    