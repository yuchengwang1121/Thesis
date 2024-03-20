import numpy as np
from pandas import read_csv

CAM = np.eye(64)
s_filename = "../MyNeuro/data/"

# write CAM
# np.savetxt(s_filename + "CAM.csv", CAM , delimiter=",",fmt='%d')

# write LUT
LUT = np.zeros((16,32))
readLUT = read_csv('../Data/LUTbin.csv',header=None).values

for i in range(len(readLUT)):
    letter = [int(x) for x in readLUT[i][0]]
    print(readLUT[i][0], "with int => ", np.array(letter))
    
    if(i>0):
        LUT = np.vstack([LUT,np.array(letter)])
    else:
        LUT = np.array(letter)
    # np.append(LUT, np.array(letter))
print(LUT)
# np.savetxt(s_filename + "LUT.csv", LUT , delimiter=",",fmt='%d')