import numpy as np
import math
from pandas import read_csv
import torch as tr
import torch.nn.functional as F
from argparse import ArgumentParser
import sys
 
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
parser.add_argument('--src1', type=str,
                        help='The path where store {src1}.csv file')
parser.add_argument('--src2', type=str,
                        help='The path where store {src2}.csv file')

args = parser.parse_args()
print("====> The parser are setting as below\n")
print("==> The path is : ", args.Path)
if not(args.STAR or args.TransSeg or args.CalculMSE):
    print("Please choose one of the mode below as 'true' ")
    print("1. STAR")
    print("2. TransSeg")
    print("3. CalculMSE")
    sys.exit(1)

print("==> The mode is : ", 'STAR' if args.STAR else 'TransSeg' if args.TransSeg else 'Calculate MSE')
if(args.CalculMSE):
    if not (args.src1 and args.src2):
        print("--src1 and --src2 are required when CalculMSE is True")
        sys.exit(1)
    else:
        print("==> The src1 is : ", args.Path+args.src1)
        print("==> The src2 is : ", args.Path+args.src2)


#######################     Read file     #######################
print("\n####################### Read file #######################\n")
q = read_csv(args.Path+'Weight0_Q.csv',header=None).values
k = read_csv(args.Path+'Weight0_K.csv',header=None).values
v = read_csv(args.Path+'Weight0_V.csv',header=None).values

In = read_csv(args.Path+'Input.csv',header=None).values
Seg_array = np.split(In, args.Seg, axis=0)
Seg_size = int(In.shape[0]/args.Seg)

print("==> The In size is ", In.shape, " with num of Seg", len(Seg_array))

#######################     The Main Function     #######################
# STAR
if(args.STAR):
    Q = np.matmul(In, q)
    K = np.matmul(In, k)
    V = np.matmul(In, v)

    scores = np.around(np.matmul(Q ,K.transpose()))
    scores = F.softmax(tr.tensor(scores), dim=-1)
    scores = np.matmul(scores, V)
    print("\n####################### STAR #######################\n")
    print("==> Saving the Softmax result of 'STAR' into the specify folder : ", args.Path, " as Soft_STAR.csv")
    np.savetxt(args.Path + "Soft_STAR.csv", scores , delimiter=",",fmt='%10.5f')

# TransSeg
elif (args.TransSeg):
    V = np.matmul(In, v)
    QK_seg = np.empty((args.Seg, Seg_size, In.shape[1]))
    for i in range(args.Seg):
        QK_seg[i] = np.matmul(np.matmul(Seg_array[i], q),k.transpose())

## Ring here
    TI = np.empty((args.Seg, In.shape[1], Seg_size))
    S  = np.empty((args.Seg, Seg_size, Seg_size))
    merge_S = []
    for i in range(args.Seg): 
        TI[i] = np.array([Seg_array[i].transpose()])

    for i in range(args.Seg):
        for j in range(args.Seg):
            result = F.softmax(tr.tensor(np.matmul(QK_seg[i],TI[j])), dim=-1)
            S[j] = np.around(result.numpy(), decimals=5)
        
        merge_S.append(np.concatenate(S, axis=1))
    Res = np.concatenate(merge_S, axis=0)
    scores = np.matmul(Res, V)

## Save result
    print("\n####################### TransSeg #######################\n")
    print("==> Saving the Softmax result of 'TransSeg' into the specify folder : ", args.Path , " as ", "Soft_TransSeg_"+ str(args.Seg) +".csv")
    np.savetxt(args.Path + "Soft_TransSeg_"+ str(args.Seg) +".csv", scores , delimiter=",",fmt='%10.5f')

#######################     Calculate MSE     #######################
if(args.CalculMSE):
    ori = read_csv(args.Path+args.src1,header=None).values
    cmp = read_csv(args.Path+args.src2,header=None).values

    mse = np.mean((ori-cmp)**2)
    print("==> The MSE result is : ", mse)

    