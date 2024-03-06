import numpy as np
import math
from pandas import read_csv
import torch as tr
import torch.nn.functional as F
 

Retransformer = False
# Retransformer = True

d_k = 128
Seg = 4
s_filename = "../Data/"

q = read_csv('../Data/Weight0_Q.csv',header=None).values
k = read_csv('../Data/Weight0_K.csv',header=None).values
v = read_csv('../Data/Weight0_V.csv',header=None).values

In = read_csv('../Data/Input.csv',header=None).values
In_1,In_2,In_3,In_4 = np.split(In, Seg, axis=0)


print("=>The In size is ", In.shape)
print("==> The q.weight size is ", q.shape)
print("==> The k.weight size is ", k.shape)
print("==> The v.weight size is ", v.shape)

if(Retransformer):
    Q = np.matmul(In, q)
    K = np.matmul(In, k)
    V = np.matmul(In, v)

    scores = np.matmul(Q ,K.transpose())
    scores = F.softmax(tr.tensor(scores), dim=-1)
    np.savetxt(s_filename + "Soft.csv", scores , delimiter=",",fmt='%10.5f')
else:
    Q_1 = np.matmul(In_1, q)
    Q_2 = np.matmul(In_2, q)
    Q_3 = np.matmul(In_3, q)
    Q_4 = np.matmul(In_4, q)

    QK_1 = np.matmul(Q_1 ,k.transpose())
    QK_2 = np.matmul(Q_2 ,k.transpose())
    QK_3 = np.matmul(Q_3 ,k.transpose())
    QK_4 = np.matmul(Q_4 ,k.transpose())

    s_1 = np.matmul(QK_1 ,In_1.transpose())
    s_2 = np.matmul(QK_2 ,In_2.transpose())
    s_3 = np.matmul(QK_3 ,In_3.transpose())
    s_4 = np.matmul(QK_4 ,In_4.transpose())

    S_1 = F.softmax(tr.tensor(s_1), dim=-1)
    S_2 = F.softmax(tr.tensor(s_2), dim=-1)
    S_3 = F.softmax(tr.tensor(s_3), dim=-1)
    S_4 = F.softmax(tr.tensor(s_4), dim=-1)

    np.savetxt(s_filename + "Soft_1.csv", S_1 , delimiter=",",fmt='%10.5f')
    np.savetxt(s_filename + "Soft_2.csv", S_2 , delimiter=",",fmt='%10.5f')
    np.savetxt(s_filename + "Soft_3.csv", S_3 , delimiter=",",fmt='%10.5f')
    np.savetxt(s_filename + "Soft_4.csv", S_4 , delimiter=",",fmt='%10.5f')