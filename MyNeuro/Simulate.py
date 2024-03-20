from subprocess import call
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--mode', action='store',  default='MyNeuro',
                        help='run different Structure in NeuroSim')
args = parser.parse_args()
if args.mode == "Retransformer" :
    call(["/bin/bash", "./data/Retransformer.sh"])								## Print the result here
elif args.mode == "TransSeg" :
    call(["/bin/bash", "./data/TransSeg.sh"])								## Print the result here
elif args.mode == "STAR" :
    call(["/bin/bash", "./data/STAR.sh"])								## Print the result here
else :
    call(["/bin/bash", "./data/trace_command.sh"])	