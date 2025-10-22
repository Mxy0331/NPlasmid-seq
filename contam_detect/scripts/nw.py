# -*- coding: utf-8 -*-
# @Time    : 2022/11/04 15:59
# @Author  : zzk
# @File    : nw.py
# Description:

# -*- coding: UTF-8 -*-
import numpy as np
#网上查清楚NW算法原理，赋分原则等等。
def dynamic_edit(str1, str2):
    len_str1 = len(str1) + 1
    len_str2 = len(str2) + 1
    newstr1 = newstr2 = ""
    matrix = np.zeros([len_str1, len_str2])
    matrix[0] = np.arange(len_str2)
    matrix[:, 0] = np.arange(len_str1)
    for i in range(1, len_str1):
        for j in range(1, len_str2):
            matrix[i][j] = min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1,
                               matrix[i - 1][j - 1] + 2*(int)(str1[i - 1] != str2[j - 1]))
    i = len_str1-1
    j = len_str2-1
    while i != 0 or j != 0:
        if (matrix[i][j] == matrix[i - 1][j - 1] and str1[i - 1] == str2[j - 1]) or matrix[i][j] == matrix[i - 1][j - 1] + 2:
            newstr1 = str1[i - 1].__add__(newstr1)
            newstr2 = str2[j - 1].__add__(newstr2)
            i = i-1
            j = j-1
        elif matrix[i][j] == matrix[i - 1][j] + 1:
            newstr1 = str1[i - 1].__add__(newstr1)
            newstr2 = "*".__add__(newstr2)
            i = i-1
        elif matrix[i][j] == matrix[i][j - 1] + 1:
            newstr1 = "*".__add__(newstr1)
            newstr2 = str2[j - 1].__add__(newstr2)
            j = j-1
    return [newstr1,newstr2]


