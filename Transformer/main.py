from torchtext.data.utils import get_tokenizer
from torchtext.vocab import build_vocab_from_iterator
from torchtext.datasets import multi30k, Multi30k
from typing import Iterable, List
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from timeit import default_timer as timer
from torch.utils.data import DataLoader
from torch.nn.utils.rnn import pad_sequence
# import from other python file
from Mask import create_masks
from Model import Transformer, get_model
from Mask import create_masks
from Translate import translate_sentence

"""
multi30k generate and token_transform
"""
# We need to modify the URLs for the dataset since the links to the original dataset are broken
# Refer to https://github.com/pytorch/text/issues/1756#issuecomment-1163664163 for more info
multi30k.URL["train"] = "https://raw.githubusercontent.com/neychev/small_DL_repo/master/datasets/Multi30k/training.tar.gz"
multi30k.URL["valid"] = "https://raw.githubusercontent.com/neychev/small_DL_repo/master/datasets/Multi30k/validation.tar.gz"

SRC_LANGUAGE = 'de'
TGT_LANGUAGE = 'en'

# Place-holders
token_transform = {}
vocab_transform = {}

token_transform[SRC_LANGUAGE] = get_tokenizer('spacy', language='de_core_news_sm')
token_transform[TGT_LANGUAGE] = get_tokenizer('spacy', language='en_core_web_sm')

# helper function to yield list of tokens
def yield_tokens(data_iter: Iterable, language: str) -> List[str]:
    language_index = {SRC_LANGUAGE: 0, TGT_LANGUAGE: 1}

    for data_sample in data_iter:
        yield token_transform[language](data_sample[language_index[language]])

# Define special symbols and indices
UNK_IDX, PAD_IDX, BOS_IDX, EOS_IDX = 0, 1, 2, 3
# Make sure the tokens are in order of their indices to properly insert them in vocab
special_symbols = ['<unk>', '<pad>', '<bos>', '<eos>']

for ln in [SRC_LANGUAGE, TGT_LANGUAGE]:
    # Training data Iterator
    train_iter = Multi30k(split='train', language_pair=(SRC_LANGUAGE, TGT_LANGUAGE))
    # Create torchtext's Vocab object
    vocab_transform[ln] = build_vocab_from_iterator(yield_tokens(train_iter, ln),
                                                    min_freq=1,
                                                    specials=special_symbols,
                                                    special_first=True)

# Set ``UNK_IDX`` as the default index. This index is returned when the token is not found.
# If not set, it throws ``RuntimeError`` when the queried token is not found in the Vocabulary.
for ln in [SRC_LANGUAGE, TGT_LANGUAGE]:
  vocab_transform[ln].set_default_index(UNK_IDX)
"""
Generate transformer model
"""
torch.manual_seed(0)
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
SRC_VOCAB_SIZE = len(vocab_transform[SRC_LANGUAGE])
TGT_VOCAB_SIZE = len(vocab_transform[TGT_LANGUAGE])
D_MODEL = 128
NHEAD = 8
FFN_HID_DIM = 128
BATCH_SIZE = 32
N_LAYERS = 3
DROPOUT = 0.1
PreTrain = None     # Pretrain if trained

model = get_model(D_MODEL, N_LAYERS, NHEAD, DROPOUT, SRC_VOCAB_SIZE, TGT_VOCAB_SIZE, PreTrain)

transformer = model.to(DEVICE)

loss_fn = torch.nn.CrossEntropyLoss(ignore_index=PAD_IDX)

optimizer = torch.optim.Adam(transformer.parameters(), lr=0.0001, betas=(0.9, 0.98), eps=1e-9)

# ``src`` and ``tgt`` language text transforms to convert raw strings into tensors indices
text_transform = {}
for ln in [SRC_LANGUAGE, TGT_LANGUAGE]:
    text_transform[ln] = sequential_transforms(token_transform[ln], #Tokenization
                                               vocab_transform[ln], #Numericalization
                                               tensor_transform) # Add BOS/EOS and create tensor

def collate_fn(batch):
    src_batch, tgt_batch = [], []
    for src_sample, tgt_sample in batch:
        src_batch.append(text_transform[SRC_LANGUAGE](src_sample.rstrip("\n")))
        tgt_batch.append(text_transform[TGT_LANGUAGE](tgt_sample.rstrip("\n")))

    src_batch = pad_sequence(src_batch, padding_value=PAD_IDX)
    tgt_batch = pad_sequence(tgt_batch, padding_value=PAD_IDX)
    return src_batch, tgt_batch

"""
Training
"""
def train_epoch(model, optimizer):
    model.train()
    losses = 0
    train_iter = Multi30k(split='train', language_pair=(SRC_LANGUAGE, TGT_LANGUAGE))
    train_dataloader = DataLoader(train_iter, batch_size=BATCH_SIZE, collate_fn=collate_fn)
    train_total = len(train_dataloader)
    
    for src, tgt in train_dataloader:
        src = src.to(DEVICE)
        tgt = tgt.to(DEVICE)
        tgt_input = tgt[:-1, :]

        src_mask, tgt_mask, src_padding_mask, tgt_padding_mask = create_masks(src, tgt_input, device=DEVICE, pad_idx=PAD_IDX)

        logits = model(src, tgt_input, src_mask, tgt_mask,src_padding_mask, tgt_padding_mask, src_padding_mask)

        optimizer.zero_grad()

        tgt_out = tgt[1:, :]
        loss = loss_fn(logits.reshape(-1, logits.shape[-1]), tgt_out.reshape(-1))
        loss.backward()
        optimizer.step()
        losses += loss.item()

    return losses / train_total

"""
Evaluate
"""
def evaluate(model):
    model.eval()
    losses = 0

    val_iter = Multi30k(split='valid', language_pair=(SRC_LANGUAGE, TGT_LANGUAGE))
    val_dataloader = DataLoader(val_iter, batch_size=BATCH_SIZE, collate_fn=collate_fn)
    val_total = len(val_dataloader)

    for src, tgt in val_dataloader:
        src = src.to(DEVICE)
        tgt = tgt.to(DEVICE)

        tgt_input = tgt[:-1, :]

        src_mask, tgt_mask, src_padding_mask, tgt_padding_mask = create_masks(src, tgt_input, device=DEVICE, pad_idx=PAD_IDX)

        logits = model(src, tgt_input, src_mask, tgt_mask,src_padding_mask, tgt_padding_mask, src_padding_mask)

        tgt_out = tgt[1:, :]
        loss = loss_fn(logits.reshape(-1, logits.shape[-1]), tgt_out.reshape(-1))
        losses += loss.item()

    return losses / val_total

"""
Run train & eval
"""
temp_loss = 999
NUM_EPOCHS = 20

for epoch in range(1, NUM_EPOCHS+1):
    start_time = timer()
    train_loss = train_epoch(transformer, optimizer)
    end_time = timer()
    val_loss = evaluate(transformer)
    if(val_loss < temp_loss):
        temp_loss = val_loss
        torch.save(model.state_dict(), 'Pretrain/model_weights.pt')
    print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Val loss: {val_loss:.3f}, "f"Epoch time = {(end_time - start_time):.3f}s"))

"""
Load module & write into .csv file
"""
# PreTrain = "Pretrain"
# WRITE_INPUT = True
# input_string = "Zwei junge weiße Männer sind im, Freien in der Nähe vieler Büsche."
# pre_transformer = get_model(D_MODEL, N_LAYERS, NHEAD, DROPOUT, SRC_VOCAB_SIZE, TGT_VOCAB_SIZE, PreTrain)

### write structure into .txt
# f = open("Model_Structure.txt", 'w')
# for p in transformer.state_dict():
#     f.write("{}\t{}\n".format(p,transformer.state_dict()[p].shape))
# f.close()

### write weight into .csv
### transformer.encoder.layers.0.self_attn.in_proj_weight => [0]
### name, weight => [1] 
### row => [0:384]
# weight_matrix = list(pre_transformer.named_parameters())[0][1].detach().numpy()
# np.savetxt("../data/weight.csv", weight_matrix, delimiter=",",fmt='%10.5f')

### split weight into Q,K,V.csv, notice to delete "header=False" in main.py/weight.csv to execute below function
# f = pd.read_csv("../data/weight.csv",header=None).values
# step = 128
# for start in range(0,len(f),step):
#     stop = start + step
#     print("The start is : ", start, " with stop ", stop)
#     if(start == 0):
#         filename = "../data/Q_weight.csv"
#     elif(start == 127):
#         filename = "../data/K_weight.csv"
#     else:
#         filename = "../data/V_weight.csv"
#     np.savetxt(filename, f[start:stop] , delimiter=",",fmt='%10.5f')

### split input into In1,In2,In3,In4.csv, notice to delete "header=False" in seq2seq.py/encode/input.csv to execute below function
# f = pd.read_csv("../data/input.csv",header=None).values
# step = 4
# for start in range(0,len(f),step):
#     stop = start + step
#     print("The start is : ", start, " with stop ", stop)
#     if(start == 0):
#         filename = "../data/In1.csv"
#     elif(start == 4):
#         filename = "../data/In2.csv"
#     elif(start == 8):
#         filename = "../data/In3.csv"
#     else:
#         filename = "../data/In4.csv"
#     np.savetxt(filename, f[start:stop] , delimiter=",",fmt='%10.5f')

### Tanspose the segment
# f = pd.read_csv("../MyNeuro/data/Input.csv",header=None).values.transpose()
# print(f.shape)
# filename = "../MyNeuro/data/TransInput.csv"
# np.savetxt(filename, f , delimiter=",",fmt='%10.5f')
# f = pd.read_csv("..../MyNeuro/data/SegIn.csv",header=None).values.transpose()
# print(f.shape)
# filename = "../MyNeuro/data/TransSegIn.csv"
# np.savetxt(filename, f , delimiter=",",fmt='%10.5f')
"""
Translate
"""
# print("The input string is : ", input_string)
# # Two young, White males are outside near many bushes.
# print("The output string is : ", translate_sentence(input_string, pre_transformer, 
#         SRC = text_transform[SRC_LANGUAGE], TGT = text_transform[TGT_LANGUAGE], devide=DEVICE))