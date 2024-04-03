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
import string
import random

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
    sentence = Variable(torch.LongTensor([indexed])).to(opt.device)

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
    parser.add_argument('-k', type=int, default=3)              # The max time find eos_tok
    parser.add_argument('-max_len', type=int, default=80)       # Max len to store stence
    parser.add_argument('-d_model', type=int, default=128)
    parser.add_argument('-n_layers', type=int, default=5)
    parser.add_argument('-src_lang', required=True)
    parser.add_argument('-trg_lang', required=True)
    parser.add_argument('-heads', type=int, default=1)
    parser.add_argument('-dropout', type=int, default=0.1)
    parser.add_argument('-no_cuda', action='store_true')
    parser.add_argument('-floyd', action='store_true')
    parser.add_argument('-sim_ratio', type=float, default=0.8)
    parser.add_argument('-test_times', type=int, default=1)
    parser.add_argument('-Writedata', action='store', default=False)
    parser.add_argument('-Segnum', type=int, default=0)
    
    opt = parser.parse_args()
    opt.device = 'cuda' if opt.no_cuda is False else 'cpu'
    if opt.device == 'cuda':
        assert torch.cuda.is_available()
    print("\n########## The parameter ##########\n")
    print("==> The device is : ", opt.device)
    print("==> The Segnum is : ", opt.Segnum)
    print("==> The test_times is : ", opt.test_times)
    print("==> The Writedata is : ", opt.Writedata)
    print("\n########## End of parameter ##########\n")

    assert opt.k > 0
    assert opt.max_len > 10

    SRC, TRG = create_fields(opt)
    model = get_model(opt, len(SRC.vocab), len(TRG.vocab))
    
    avg_acc = 0

    with open('./data/english.txt', 'r') as EF:
        English_content = EF.readlines()
    with open('./data/french.txt', 'r') as FF:
        French_content = FF.readlines()

    print("The len of data is =>",len(English_content))

    
    # opt.text =  "Between meals, he usually manages to stow away a generous supply of candy, ice cream, popcorn and fruit."
    # Ans = "Entre les repas, il s'arrange d'ordinaire pour mettre de côté une abondante réserve de sucreries, de crème glacée, de pop-corn et de fruits."
    
    # opt.text =  "One of the best ways to help us is to translate from a foreign language you know into your own native language or strongest language."
    # Ans = "L'une des meilleures manières de nous aider est de traduire d'une langue étrangère que vous connaissez vers votre propre langue natale, ou la plus forte de vos langues."
    # phrase = translate(opt, model, SRC, TRG)
    # print("After translate :", phrase)
    # print("After Answer :", Ans)
    # res = calculate_equal_ratio(preprocess_string(phrase), preprocess_string(Ans.strip()))
    # print("HERE", res)
    
    # Calculate the accuracy
    for t in range(opt.test_times):
        random_numbers = random.sample(range(len(English_content)), 100)
        hit = 0
        for i in range(len(random_numbers)):
            opt.text =  English_content[random_numbers[i]].strip()
            phrase = translate(opt, model, SRC, TRG)
            # print("After translate :", phrase)
            # print("After Answer :", French_content[random_numbers[i]].strip())
            if(calculate_equal_ratio(preprocess_string(phrase), preprocess_string(French_content[random_numbers[i]].strip())) >= opt.sim_ratio):
                hit += 1
        acc = round((hit/len(random_numbers))*100, 3)
        print("The accuracy in %d times is %.2f" %(t, acc))
        avg_acc += acc
    print("The overall average accuracy is %.2f%%" %round((avg_acc/opt.test_times), 3))

                
def preprocess_string(s):
    s = s.translate(str.maketrans('', '', string.punctuation))
    return s

def calculate_equal_ratio(string1, string2):
    count = 0
    min_length = min(len(string1.split()), len(string2.split()))
    equal_min_length = min(len(string1), len(string2))
    ans_split = string2.split()
    for i in range(len(ans_split)):
        index = string1.find(ans_split[i])
        # print("The ans is ", ans_split[i], "w i ", index)
        if(index != -1): count = count + 1

    equal_count = sum(1 for char1, char2 in zip(string1, string2) if char1 == char2)

    equal_ratio = max(count/min_length, equal_count/equal_min_length)  if min_length > 0 else 0
    # print("First ==> ", count/min_length, " with count : ", count , " & length : ", min_length)
    # print("Second ==> ", equal_count/equal_min_length, " with equal_count : ", equal_count , " & length : ", equal_min_length)
    # print(equal_ratio)
    return equal_ratio

if __name__ == '__main__':
    main()
