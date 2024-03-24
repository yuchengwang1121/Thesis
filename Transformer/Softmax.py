import numpy as np
import math
from pandas import read_csv
import torch as tr
import torch.nn.functional as F
from argparse import ArgumentParser
 
#######################     The Parser     #######################
parser = ArgumentParser()
parser.add_argument('--Path', default="../Data/",
                        help='The path to store .csv file')
parser.add_argument('--STAR', action='store', default=False,
                        help='Run the Softmax on STAR')
parser.add_argument('--TransSeg', action='store', default=False,
                        help='Run the Softmax on TransSeg')
parser.add_argument('--Seg', type=int,  default=4,
                        help='The number of segment')
parser.add_argument('--CalculMSE', action='store', default=False,
                        help='Calculate MSE between different number of segment')

args = parser.parse_args()
s_filename = args.Path

#######################     Read file     #######################
q = read_csv(args.Path+'Weight0_Q.csv',header=None).values
k = read_csv(args.Path+'Weight0_K.csv',header=None).values
v = read_csv(args.Path+'Weight0_V.csv',header=None).values

In = read_csv(args.Path+'Input.csv',header=None).values
Seg_array = np.split(In, args.Seg, axis=0)
print("=>The In size is ", In.shape, " with num of Seg", len(Seg_array))
print("==> The q.weight size is ", q.shape)
print("==> The k.weight size is ", k.shape)
print("==> The v.weight size is ", v.shape)

#######################     The Main Function     #######################
# STAR
if(args.STAR):
    Q = np.matmul(In, q)
    K = np.matmul(In, k)
    V = np.matmul(In, v)

    scores = np.around(np.matmul(Q ,K.transpose()))
    scores = F.softmax(tr.tensor(scores), dim=-1)
    print("==> Saving the Softmax result of 'STAR' into the specify folder : ", s_filename)
    np.savetxt(s_filename + "Soft_STAR.csv", scores , delimiter=",",fmt='%10.5f')

# TransSeg
elif (args.TransSeg):
    QK_seg = np.empty((args.Seg, args.Seg, In.shape[1]))
    for i in range(args.Seg):
        QK_seg[i] = np.matmul(np.matmul(Seg_array[i], q),k.transpose())

## Ring here
    TI = np.empty((args.Seg, In.shape[1], args.Seg))
    S  = np.empty((args.Seg, args.Seg, args.Seg))
    merge_S = []
    for i in range(args.Seg): 
        TI[i] = np.array([Seg_array[i].transpose()])

    for i in range(args.Seg):
        for j in range(args.Seg):
            result = F.softmax(tr.tensor(np.matmul(QK_seg[i],TI[j])), dim=-1)
            S[j] = np.around(result.numpy(), decimals=5)
        
        merge_S.append(np.concatenate(S, axis=1))
    Res = np.concatenate(merge_S, axis=0)

## Save result
    print("==> Saving the Softmax result of 'TransSeg' into the specify folder : ", s_filename)
    np.savetxt(s_filename + "Soft__TransSeg.csv", Res , delimiter=",",fmt='%10.5f')