from deep_translator import GoogleTranslator
from argparse import ArgumentParser
import time

parser = ArgumentParser()
parser.add_argument('--pivot', type=int,  required=True,
                        help='The pivot to start')
parser.add_argument('--TGT_lan', default="None",
                        help='The path to store .csv file')

args = parser.parse_args()
src_path = "./data/english.txt"
dst_path = "./data/"+ str(args.TGT_lan) + ".txt"
steps = 10
start = 10000*args.pivot
end   = start + (steps*1000)
start_time = time.time()


with open(src_path, "r", encoding="utf-8") as file, \
     open(dst_path, "a", encoding="utf-8") as o_file:

    for index, line in enumerate(file):
        if(start <= index and index < end):
            translated = GoogleTranslator(source='auto', target='de').translate(line.strip())
            o_file.write(translated + "\n")
            if(index % 1000 == 0):
                print("{}m : Now traslate from {} to {}".format((time.time() - start_time)//60,index, index+1000))
        elif(index > end):
            break
