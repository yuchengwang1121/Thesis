import numpy as np
import math
from pandas import read_csv
import torch as tr
import torch.nn.functional as F
from argparse import ArgumentParser
import sys
 
#######################     The Parser     #######################
parser = ArgumentParser()
# Path 
parser.add_argument('--WeightPath', default="./data/Weight/",
                        help='The path to store .csv file')
parser.add_argument('--InputPath', default="./data/Input/",
                        help='The path to store .csv file')
parser.add_argument('--ResultPath', default="./data/Result/",
                        help='The path to store .csv file')
# Generater Attention score
parser.add_argument('--Layer', type=int,  default=None,
                        help='The Layer to simulate.')
parser.add_argument('--STAR', action='store', default=False,
                        help='Run the Softmax on STAR')
parser.add_argument('--TransSeg', action='store', default=False,
                        help='Run the Softmax on TransSeg')
parser.add_argument('--Seg', type=int,  default=4,
                        help='The number of segment')
# Calculate Mean Square Error
parser.add_argument('--CalculMSE', action='store', default=False,
                        help='Calculate MSE between different number of segment')
parser.add_argument('--src1', type=str,default=None,
                        help='The path where store {src1}.csv file')
parser.add_argument('--src2', type=str,default=None,
                        help='The path where store {src2}.csv file')
# Generate CAM&LUT data for NeuroSim
parser.add_argument('--GenCAMLUT', action='store', default=False,
                        help='Calculate MSE between different number of segment')
parser.add_argument('--CAM_size', type=int,  default=64,
                        help='The Size of Match Vector')
parser.add_argument('--write_CAM', action='store', default=False,
                        help='Write the CAM into .csv file')
parser.add_argument('--LUT_size', type=int,  default=16,
                        help='The Size of Match Vector')
parser.add_argument('--write_LUT', action='store', default=False,
                        help='Write the LUT to fit NeuroSim')
parser.add_argument('--Bin_Pre',  type=int,  default=32,
                        help='The Precise to store LUT into binary form')
# Generate Input data for NeuroSim
parser.add_argument('--TransSeg_Arch', action='store', default=False,
                        help='Calculate MSE between different number of segment')

args = parser.parse_args()
print("====> The parser are setting as below\n")
print("#######################  Path setting  #######################")
print("==> The weight path is : ", args.WeightPath)
print("==> The input  path is : ", args.InputPath)
print("==> The result path is : ", args.ResultPath)

if not(args.STAR or args.TransSeg or args.CalculMSE or args.GenCAMLUT or args.TransSeg_Arch):
    print("Please choose one of the mode below as 'true'.")
    print("1. STAR")
    print("2. TransSeg")
    print("3. CalculMSE")
    print("4. GenCAMLUT")
    print("5. TransSeg_Arch")
    sys.exit(1)

print("\n#######################  Mode setting  #######################")
print("==> The mode is : ", 'STAR' if args.STAR else 'TransSeg' if args.TransSeg else 'Calculate MSE'
                                   if args.CalculMSE else 'Generate CAM&LUT' if args.GenCAMLUT else 'Generate Input')

if(args.TransSeg_Arch):
    print("==> The Segment number is : ", args.Seg)
elif(args.GenCAMLUT):
    print("==> The CAM size is : ", args.CAM_size)
    print("==> The LUT size is : ", args.LUT_size, " with precision ", args.Bin_Pre)
elif(args.CalculMSE):
    if (args.src1==None or  args.src2==None):
        print("--src1 and --src2 are required when CalculMSE is True.")
        sys.exit(1)
    else:
        print("==> The src1 is : ", args.ResultPath+args.src1)
        print("==> The src2 is : ", args.ResultPath+args.src2)
else:
    if (args.Layer == None):
        print("Please choose the Layer.")
        sys.exit(1)
    else:
        # Read file
        print("\n####################### Read file #######################")
        q = read_csv(args.WeightPath+'weight' + str(args.Layer) + '_Q.csv',header=None).values
        k = read_csv(args.WeightPath+'weight' + str(args.Layer) + '_K.csv',header=None).values
        v = read_csv(args.WeightPath+'weight' + str(args.Layer) + '_V.csv',header=None).values

        In = read_csv(args.InputPath+'input_' + str(args.Layer) + '.csv',header=None).values
        Seg_array = np.split(In, args.Seg, axis=0)
        Seg_size = int(In.shape[0]/args.Seg)


#######################     The Sub Function     #######################
# Change LUT2BIN
def float_to_binary(num, precision=args.Bin_Pre):
    integer_part = int(num)
    fractional_part = num - integer_part

    integer_binary = bin(integer_part)[2:]
    fractional_binary = ""

    for _ in range(precision):
        fractional_part *= 2
        bit = int(fractional_part)
        fractional_binary += str(bit)
        fractional_part -= bit

    binary_representation = integer_binary + fractional_binary
    return binary_representation

#######################     The Main Function     #######################
# STAR
if(args.STAR):
    Q = np.matmul(In, q)
    K = np.matmul(In, k)
    V = np.matmul(In, v)

    scores = np.around(np.matmul(Q ,K.transpose()))
    scores = F.softmax(tr.tensor(scores), dim=-1)
    attscore = np.matmul(scores, V)
    print("\n####################### Saving STAR #######################\n")
    print("==> Saving the Softmax result of 'STAR' into the specify folder : ", args.ResultPath, "Layer" , str(args.Layer), " as Attscore_STAR.csv")
    np.savetxt(args.ResultPath + "Layer" + str(args.Layer) + "/Attscore_STAR.csv", attscore , delimiter=",",fmt='%10.5f')
    print("==> Saving the Softmax result of 'STAR' into the specify folder : ", args.ResultPath, "Layer" , str(args.Layer), " as Soft_STAR.csv")
    np.savetxt(args.ResultPath + "Layer" + str(args.Layer) + "/Soft_STAR.csv", scores , delimiter=",",fmt='%10.5f')

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
    attscore = np.matmul(Res, V)

## Save result
    print("\n####################### Saving TransSeg #######################\n")
    print("==> Saving the Softmax result of 'TransSeg' into the specify folder : ", args.ResultPath , " as ", "Attscore_TransSeg_"+ str(args.Seg) +".csv")
    np.savetxt(args.ResultPath + "Layer" + str(args.Layer) + "/Attscore_TransSeg_" + str(args.Seg) +".csv", attscore , delimiter=",",fmt='%10.5f')
    print("==> Saving the Softmax result of 'TransSeg' into the specify folder : ", args.ResultPath , " as ", "Soft_TransSeg_"+ str(args.Seg) +".csv")
    np.savetxt(args.ResultPath + "Layer" + str(args.Layer) + "/Soft_TransSeg_" + str(args.Seg) +".csv", Res , delimiter=",",fmt='%10.5f')

#######################       Calculate MSE       #######################
if(args.CalculMSE):
    ori = read_csv(args.ResultPath+args.src1,header=None).values
    cmp = read_csv(args.ResultPath+args.src2,header=None).values

    mse = np.mean((ori-cmp)**2)
    print("==> The MSE result is : ", mse)

#######################          Gen LUT          #######################
if(args.GenCAMLUT):
    # write CAM
    print("\n####################### Gen CAM #######################\n")
    s_filename = './data/CAM_LUT/'
    CAM = np.eye(args.CAM_size)
    if(args.write_CAM):
        print("==> Saving the CAM into the specify folder : ", s_filename)
        np.savetxt(s_filename + "CAM.csv", CAM , delimiter=",",fmt='%d')

    # write LUT to Binary form
    print("\n####################### Transfer from FP to Bin #######################\n")
    lut = np.zeros(args.LUT_size)
    for i in range(-args.LUT_size+1,1,1):
        lut[i+args.LUT_size-1] = format(math.exp(i), '.8f')  #get precision 8 
    with open(s_filename + "LUTbin.csv", "w") as f:
        for i in range(len(lut)):
            s = float_to_binary(lut[i])
            f.write(s + "\n")

    # Split the LUTbin to fit NeuroSim
    LUT = np.zeros((args.LUT_size,args.Bin_Pre))
    readLUT = read_csv(s_filename + 'LUTbin.csv',header=None).values

    for i in range(len(readLUT)):
        letter = [int(x) for x in readLUT[i][0]]
        if(i>0):
            LUT = np.vstack([LUT,np.array(letter)])
        else:
            LUT = np.array(letter)
    print("\n####################### Gen LUT #######################\n")
    if(args.write_LUT):
        print("==> Saving the LUT into the specify folder : ", s_filename)
        np.savetxt(s_filename + "LUT.csv", LUT , delimiter=",",fmt='%d')

#######################       Gen Data for Neuro       #######################
if(args.TransSeg_Arch):
    input = read_csv(args.InputPath+"input_0.csv",header=None).values
    segments = np.array_split(input, args.Seg, axis=0)[0]
    np.savetxt(args.InputPath + "TransSegIn_" + str(args.Seg) +".csv", segments , delimiter=",",fmt='%10.5f')
    np.savetxt(args.InputPath + "TransSegIn_T" + str(args.Seg) +".csv", segments.transpose() , delimiter=",",fmt='%10.5f')
