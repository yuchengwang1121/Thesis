from subprocess import call
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--mode', action='store',  default='MyNeuro',
                        help='run different Structure in NeuroSim')
parser.add_argument('--segnum', type=int,  default=None,
                        help='Num of Segment.')
args = parser.parse_args()
if args.mode == "TransSeg_Arch" :
    if (args.segnum == None):
        print("The segnumber is needed.")
    else:
        call(["/bin/bash", "./data/Shell/TransSeg_Arch_" + str(args.segnum) +".sh"])								## Print the result here
elif args.mode == "TransSeg_Soft" :
    call(["/bin/bash", "./data/Shell/TransSeg_Soft.sh"])								## Print the result here
elif args.mode == "STAR_Soft" :
    call(["/bin/bash", "./data/Shell/STAR_Soft.sh"])								## Print the result here
else :
    call(["/bin/bash", "./data/Shell/STAR_Arch.sh"])	