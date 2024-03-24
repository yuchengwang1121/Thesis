import numpy as np
import math
from pandas import read_csv
from argparse import ArgumentParser

#######################     The Parser     #######################
parser = ArgumentParser()
parser.add_argument('--Path', default="../MyNeuro/data/",
                        help='The path to store .csv file')
parser.add_argument('--CAM_size', type=int,  default=64,
                        help='The Size of Match Vector')
parser.add_argument('--write_CAM', action='store', default=False,
                        help='Write the CAM into .csv file')
parser.add_argument('--LUT_len', type=int,  default=16,
                        help='The Size of Match Vector')
parser.add_argument('--write_LUT', action='store', default=False,
                        help='Write the LUT to fit NeuroSim')
parser.add_argument('--Bin_Pre',  type=int,  default=32,
                        help='The Precise to store LUT into binary form')

args = parser.parse_args()
s_filename = args.Path
print("The File will be store into ", s_filename)

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
# write CAM
CAM = np.eye(args.CAM_size)
if(args.write_CAM):
    print("==> Saving the CAM into the specify folder : ", s_filename)
    np.savetxt(s_filename + "CAM.csv", CAM , delimiter=",",fmt='%d')

# write LUT to Binary form
lut = np.zeros(args.LUT_len)
for i in range(-args.LUT_len+1,1,1):
    lut[i+args.LUT_len-1] = format(math.exp(i), '.8f')  #get precision 8 
with open("../Data/LUTbin.csv", "w") as f:
    for i in range(len(lut)):
        s = float_to_binary(lut[i])
        f.write(s + "\n")


# Split the LUTbin to fit NeuroSim
LUT = np.zeros((args.LUT_len,args.Bin_Pre))
readLUT = read_csv('../Data/LUTbin.csv',header=None).values

for i in range(len(readLUT)):
    letter = [int(x) for x in readLUT[i][0]]
    if(i>0):
        LUT = np.vstack([LUT,np.array(letter)])
    else:
        LUT = np.array(letter)
    if(args.write_LUT):
        print("==> Saving the LUT into the specify folder : ", s_filename)
        np.savetxt(s_filename + "LUT.csv", LUT , delimiter=",",fmt='%d')





