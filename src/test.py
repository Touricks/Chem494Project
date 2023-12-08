

import numpy as np
def levenshtein_distance(str_a,str_b):
    str_a=str_a.lower()
    str_b=str_b.lower()
    matrix_ed=np.zeros((len(str_a)+1,len(str_b)+1),dtype=int)
    matrix_ed[0]=np.arange(len(str_b)+1)
    matrix_ed[:,0] = np.arange(len(str_a) + 1)
    for i in range(1,len(str_a)+1):
        for j in range(1,len(str_b)+1):
            # 表示删除a_i
            dist_1 = matrix_ed[i - 1, j] + 1
            # 表示插入b_i
            dist_2 = matrix_ed[i, j - 1] + 1
            # 表示替换b_i
            dist_3 = matrix_ed[i - 1, j - 1] + (1 if str_a[i - 1] != str_b[j - 1] else 0)
            #取最小距离
            matrix_ed[i,j]=np.min([dist_1, dist_2, dist_3])
    print(matrix_ed[-1][-1])
    return matrix_ed[-1,-1]

def similarity2(x,y):
    str1=''.join(x)
    str2=''.join(y)
    print(str1,str2)
    dis=levenshtein_distance(str1,str2)
    if dis/max(len(str1),len(str2))<=0.2:
        return True
    else:
        return False

x=["kitten"]
y=["sitting"]
print(similarity2(x,y))