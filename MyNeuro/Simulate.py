from subprocess import call
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--mode', action='store', default='MyNeuro',
                        help='run different Structure in NeuroSim')
parser.add_argument('--segnum', type=int,  default=None,
                        help='Num of Segment.')
args = parser.parse_args()
if args.mode == "TransSeg_Arch" :
    if (args.segnum == None):
        print("The segnumber is needed.")
    else:
        with open("./data/Shell/Arch.sh", 'r') as file:
            lines = file.read().split()
        lines[0] = "./Simulator/TransSeg"
        lines[2] = str(args.segnum)
        lines[6] = "./data/TransSeg_Arch/Weight.csv"
        lines[7] = "./data/TransSeg_Arch/TransSegIn_T" + str(args.segnum) +".csv"
        lines[8] = "./data/TransSeg_Arch/Input.csv"
        lines[9] = "./data/TransSeg_Arch/TransSegIn_" + str(args.segnum) +".csv"
        lines[10] = "./data/TransSeg_Arch/Soft_TransSeg_" + str(args.segnum) +".csv"
        with open("./data/Shell/Arch.sh", 'w') as file:
            file.write(' '.join(lines))
        call(["/bin/bash", "./data/Shell/Arch.sh"])
elif args.mode == "TransSeg_Soft" :
    if (args.segnum == None):
        print("The segnumber is needed.")
    else:
        ## Modify the shell
        with open("./data/Shell/Softmax.sh", 'r') as file:
            lines = file.read().split()
        lines[2] = "0"
        lines[6] = str(args.segnum)
        lines[7] = "./data/TransSeg_Soft/CAM_Weight.csv"
        lines[8] = "./data/TransSeg_Soft/LUT_Weight_TransSeg.csv"
        lines[9] = "./data/TransSeg_Soft/CAMInput_TransSeg_" + str(args.segnum) +".csv"
        lines[9] = "./data/TransSeg_Soft/LUTInput_TransSeg.csv"
        with open("./data/Shell/Softmax.sh", 'w') as file:
            file.write(' '.join(lines))
        call(["/bin/bash", "./data/Shell/Softmax.sh"])
elif args.mode == "STAR_Soft" :
    with open("./data/Shell/Softmax.sh", 'r') as file:
        lines = file.read().split()
    lines[2] = "1"
    lines[6] = "1"
    lines[7] = "./data/STAR_Soft/CAM_Weight.csv"
    lines[8] = "./data/STAR_Soft/LUT_Weight_STAR.csv"
    lines[9] = "./data/STAR_Soft/CAMInput_STAR.csv"
    lines[10] = "./data/STAR_Soft/LUTInput_STAR.csv"
    with open("./data/Shell/Softmax.sh", 'w') as file:
        file.write(' '.join(lines))
    call(["/bin/bash", "./data/Shell/Softmax.sh"])
else :
    with open("./data/Shell/Arch.sh", 'r') as file:
        lines = file.read().split()
    lines[0] = "./Simulator/STAR"
    lines[6] = "./data/STAR_Arch/Weight.csv"
    lines[7] = "./data/STAR_Arch/TransInput.csv"
    lines[8] = "./data/STAR_Arch/Input.csv"
    lines[9] = "./data/STAR_Arch/Input.csv"
    lines[10] = "./data/STAR_Arch/Soft_STAR.csv"
    with open("./data/Shell/Arch.sh", 'w') as file:
        file.write(' '.join(lines))
    print(lines)
    call(["/bin/bash", "./data/Shell/Arch.sh"])