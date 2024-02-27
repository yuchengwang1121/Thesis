from subprocess import call
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--mode', action='store',  default='MyNeuro',
                        help='run different Structure in NeuroSim')
args = parser.parse_args()
if args.mode == "Retransformer" :
    call(["/bin/bash", "./data/Retransformer.sh"])								## Print the result here
elif args.mode == "Softmax" :
    call(["/bin/bash", "./data/softmax.sh"])								## Print the result here
else :
    call(["/bin/bash", "./data/trace_command.sh"])	