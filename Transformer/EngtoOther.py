from deep_translator import GoogleTranslator
from argparse import ArgumentParser
import time

parser = ArgumentParser()
parser.add_argument('--pivot', type=float,  required=True,
                        help='The pivot to start')
parser.add_argument('--TGT_lan',required=True,
                        help='The path to store .csv file')

args = parser.parse_args()
src_path = "./data/english.txt"
dst_path = "./data/"+ str(args.TGT_lan) + ".txt"
steps = 10
start = 10000*args.pivot
end   = start + (steps*1000)
start_time = time.time()

if(args.TGT_lan == 'german'):
    TGT = 'de'
elif(args.TGT_lan == 'spanish'):
    TGT = 'es'
elif(args.TGT_lan == 'portuguese'):
    TGT = 'pt'
elif(args.TGT_lan == 'dutch'):
    TGT = 'nl'
elif(args.TGT_lan == 'italian'):
    TGT = 'it'
else:
    print("The language is not supported")
    quit()

with open(src_path, "r", encoding="utf-8") as file, \
     open(dst_path, "a", encoding="utf-8") as o_file:

    for index, line in enumerate(file):
        if(start <= index and index < end):
            translated = GoogleTranslator(source='auto', target=TGT).translate(line.strip())
            o_file.write(translated + "\n")
            if(index % 1000 == 0):
                print("{}m : Now traslate from {} to {}".format((time.time() - start_time)//60,index, index+1000))
        elif(index > end):
            break
