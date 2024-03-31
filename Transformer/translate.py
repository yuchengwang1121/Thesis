import argparse
import time
import torch
from Function.Models import get_model
from Function.Process import *
import torch.nn.functional as F
from Function.Optim import CosineWithRestarts
from Function.Batch import create_masks
import pdb
import dill as pickle
import argparse
from Function.Models import get_model
from Function.Beam import beam_search
from nltk.corpus import wordnet
from torch.autograd import Variable
import re
# import string

def get_synonym(word, SRC):
    syns = wordnet.synsets(word)
    for s in syns:
        for l in s.lemmas():
            if SRC.vocab.stoi[l.name()] != 0:
                return SRC.vocab.stoi[l.name()]
            
    return 0

def multiple_replace(dict, text):
  # Create a regular expression  from the dictionary keys
  regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

  # For each match, look-up corresponding value in dictionary
  return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text) 

def translate_sentence(sentence, model, opt, SRC, TRG):
    
    model.eval()
    indexed = []
    sentence = SRC.preprocess(sentence)
    for tok in sentence:
        if SRC.vocab.stoi[tok] != 0 or opt.floyd is True:
            indexed.append(SRC.vocab.stoi[tok])
        else:
            indexed.append(get_synonym(tok, SRC))
    sentence = Variable(torch.LongTensor([indexed]))
    # print("==>The sentence is : ", sentence, " with size ", sentence.size())

    sentence = beam_search(sentence, model, SRC, TRG, opt)

    return multiple_replace({' ?': '?', ' !': '!', ' .': '.', '\' ': '\'', ' ,': ','}, sentence)

def translate(opt, model, SRC, TRG):
    sentences = opt.text.lower().split('.')
    translated = []
    for sentence in sentences:
        translated.append(translate_sentence(sentence + "." , model, opt, SRC, TRG).capitalize())
        break

    return (' '.join(translated))


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-load_weights', required=True)
    parser.add_argument('-k', type=int, default=3)
    parser.add_argument('-max_len', type=int, default=80)
    parser.add_argument('-d_model', type=int, default=128)
    parser.add_argument('-n_layers', type=int, default=3)
    parser.add_argument('-src_lang', required=True)
    parser.add_argument('-trg_lang', required=True)
    parser.add_argument('-heads', type=int, default=1)
    parser.add_argument('-dropout', type=int, default=0.1)
    parser.add_argument('-no_cuda', action='store_true')
    parser.add_argument('-floyd', action='store_true')
    parser.add_argument('-Segnum', type=int, default=0)
    
    opt = parser.parse_args()
    opt.device = 'cpu'
    print("The opt is", opt)
    # opt.device = 'cuda' if opt.no_cuda is False else 'cpu'
 
    assert opt.k > 0
    assert opt.max_len > 10

    SRC, TRG = create_fields(opt)
    model = get_model(opt, len(SRC.vocab), len(TRG.vocab))
    ans_pointer = 0
    hit = 0
    
    # while True:
    #     opt.text =input("Enter a sentence to translate (type 'f' to load from file, or 'q' to quit):\n")
    #     if opt.text=="q":
    #         break
    #     if opt.text=='f':
    #         fpath =input("Enter a sentence to translate (type 'f' to load from file, or 'q' to quit):\n")
    #         try:
    #             opt.text = ' '.join(open(opt.text, encoding='utf-8').read().split('\n'))
    #         except:
    #             print("error opening or reading text file")
    #             continue
    #     # print("=> Input String is : " + opt.text)
    #     phrase = translate(opt, model, SRC, TRG)
    #     print('> '+ phrase + '\n')

    

    with open('./data/english.txt', 'r') as EF:
        English_content = EF.readlines()
    with open('./data/french.txt', 'r') as FF:
        French_content = FF.readlines()

 
    for src in English_content[:9]:
        opt.text = src
        phrase = translate(opt, model, SRC, TRG)
        print("After pre=>", phrase , " with size :", len(phrase))
        print("After Ans pre=>", French_content[ans_pointer] , " with size :", len(French_content[ans_pointer]))
        if(phrase == French_content[ans_pointer]):
            hit += 1
        ans_pointer += 1

    # hex_list = [hex(ord(char)) for char in preprocess_string(phrase)]
    # hex_list2 = [hex(ord(char)) for char in preprocess_string(French_content[ans_pointer])]
    # for hex_value in hex_list:
    #     print(hex_value )
    # for hex_value in hex_list2:
    #     print(hex_value )
    # print(len(preprocess_string(phrase)))
    # print(len(preprocess_string(French_content[ans_pointer])))
    print("The hit is", hit, " with ", ans_pointer)
    print("The over accuracy is ", round((hit/ans_pointer)*100, 3))
                
# def preprocess_string(s):
#     s = s.translate(str.maketrans('', '', string.punctuation + ' '))
#     return s

if __name__ == '__main__':
    main()
