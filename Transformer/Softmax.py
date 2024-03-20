import numpy as np
import math
from pandas import read_csv
import torch as tr
import torch.nn.functional as F
 

Retransformer = False
# Retransformer = True
# LUT = True
LUT = False

d_k = 128
Seg = 4
lutlen = 15
s_filename = "../Data/"

q = read_csv('../Data/Weight0_Q.csv',header=None).values
k = read_csv('../Data/Weight0_K.csv',header=None).values
v = read_csv('../Data/Weight0_V.csv',header=None).values

In = read_csv('../Data/Input.csv',header=None).values
Qk_11 = read_csv('../Data/SQK_11.csv',header=None).values
Qk_21 = read_csv('../Data/SQK_21.csv',header=None).values
Qk_31 = read_csv('../Data/SQK_31.csv',header=None).values
Qk_41 = read_csv('../Data/SQK_41.csv',header=None).values
In_1,In_2,In_3,In_4 = np.split(In, Seg, axis=0)


print("=>The In size is ", In.shape)
print("==> The q.weight size is ", q.shape)
print("==> The k.weight size is ", k.shape)
print("==> The v.weight size is ", v.shape)

# scores = np.around(Qk_21).reshape(-1,1)
# np.savetxt(s_filename + "QK21_onecol.csv", scores ,fmt='%d')
# scores = np.around(Qk_31).reshape(-1,1)
# np.savetxt(s_filename + "QK31_onecol.csv", scores ,fmt='%d')
# scores = np.around(Qk_41).reshape(-1,1)
# np.savetxt(s_filename + "QK41_onecol.csv", scores ,fmt='%d')

# if(Retransformer):
#     Q = np.matmul(In, q)
#     K = np.matmul(In, k)
#     V = np.matmul(In, v)

#     scores = np.around(np.matmul(Q ,K.transpose())).reshape(-1,1)
#     # np.savetxt(s_filename + "QK_onecol.csv", scores ,fmt='%d')
#     scores = F.softmax(tr.tensor(scores), dim=-1)
#     # np.savetxt(s_filename + "Soft.csv", scores , delimiter=",",fmt='%10.5f')
# else:
#     Q_1 = np.matmul(In_1, q)
#     Q_2 = np.matmul(In_2, q)
#     Q_3 = np.matmul(In_3, q)
#     Q_4 = np.matmul(In_4, q)

#     QK_1 = np.matmul(Q_1 ,k.transpose())
#     QK_2 = np.matmul(Q_2 ,k.transpose())
#     QK_3 = np.matmul(Q_3 ,k.transpose())
#     QK_4 = np.matmul(Q_4 ,k.transpose())

#     ## Ring here
#     TI = np.array([In_1.transpose(), In_2.transpose(), In_3.transpose(), In_4.transpose()])
#     for i in range(len(TI)):
#         s_1 = np.matmul(QK_1 ,TI[i])
#         s_2 = np.matmul(QK_2 ,TI[(i+1)%4])
#         s_3 = np.matmul(QK_3 ,TI[(i+2)%4])
#         s_4 = np.matmul(QK_4 ,TI[(i+3)%4])

#         # np.savetxt(s_filename + "SQK_"+str(4 if((i+1)%4 == 0) else  (i+1)%4)+"1.csv", s_1 , delimiter=",",fmt='%10.5f')

#         S_1 = F.softmax(tr.tensor(s_1), dim=-1)
#         S_2 = F.softmax(tr.tensor(s_2), dim=-1)
#         S_3 = F.softmax(tr.tensor(s_3), dim=-1)
#         S_4 = F.softmax(tr.tensor(s_4), dim=-1)

#         ## Soft_(col,row).csv
#         # np.savetxt(s_filename + "Soft_"+str(4 if((i+1)%4 == 0) else  (i+1)%4)+"1.csv", S_1 , delimiter=",",fmt='%10.5f')
#         # np.savetxt(s_filename + "Soft_"+str(4 if((i+2)%4 == 0) else  (i+2)%4)+"2.csv", S_2 , delimiter=",",fmt='%10.5f')
#         # np.savetxt(s_filename + "Soft_"+str(4 if((i+3)%4 == 0) else  (i+3)%4)+"3.csv", S_3 , delimiter=",",fmt='%10.5f')
#         # np.savetxt(s_filename + "Soft_"+str(4 if((i+4)%4 == 0) else  (i+4)%4)+"4.csv", S_4 , delimiter=",",fmt='%10.5f')

def float_to_binary(num, precision=31):
    integer_part = int(num)  # 提取整数部分
    fractional_part = num - integer_part  # 提取小数部分

    integer_binary = bin(integer_part)[2:]  # 将整数部分转换为二进制
    fractional_binary = ""

    # 计算小数部分的二进制表示
    for _ in range(precision):
        fractional_part *= 2
        bit = int(fractional_part)
        fractional_binary += str(bit)
        fractional_part -= bit

    binary_representation = integer_binary + fractional_binary
    return binary_representation


# if(LUT):
#     lut = np.zeros(lutlen+1)
#     for i in range(-lutlen,1,1):
#         lut[i+lutlen] = math.exp(i)
#     print(lut)
#     np.savetxt(s_filename + "LUT.csv", lut , delimiter=",",fmt='%10.8f')
# else:
#     l = read_csv('../Data/LUT.csv',header=None).values
#     with open("../Data/LUTbin.csv", "w") as f:
#         for i in range(len(l)):
#             print("The float is ", l[i])
#             s = float_to_binary(l[i])
#             print("The bin is ", s)
#             f.write(s + "\n")
