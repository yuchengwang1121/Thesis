import torch
import torch.nn as nn
import numpy as np
from Function.Sublayers import FeedForward, MultiHeadAttention, Norm

class EncoderLayer(nn.Module):
    def __init__(self, d_model, heads, dropout=0.1):
        super().__init__()
        self.norm_1 = Norm(d_model)
        self.norm_2 = Norm(d_model)
        self.attn = MultiHeadAttention(heads, d_model, dropout=dropout)
        self.ff = FeedForward(d_model, dropout=dropout)
        self.dropout_1 = nn.Dropout(dropout)
        self.dropout_2 = nn.Dropout(dropout)
        
    def forward(self, x, mask, Writedata, Layer):
        i_filename = "./Input"
        x2 = self.norm_1(x)
        print("Writing the Norm_in input in ", i_filename + str(Layer) + ".csv", " with size ", x2.size()[1], " * ", x2.size()[2])
        np.savetxt(i_filename+ "/input_" + str(Layer) + ".csv" , x2.detach().numpy().reshape(x2.size()[1],x2.size()[2]) , delimiter=",",fmt='%10.5f')
        print("-------------------EncoderLayer %d-------------------" % Layer)
        x = x + self.dropout_1(self.attn(x2,x2,x2,mask, Writedata, Layer))
        x2 = self.norm_2(x)
        x = x + self.dropout_2(self.ff(x2))
        return x
    
# build a decoder layer with two multi-head attention layers and
# one feed-forward layer
class DecoderLayer(nn.Module):
    def __init__(self, d_model, heads, dropout=0.1):
        super().__init__()
        self.norm_1 = Norm(d_model)
        self.norm_2 = Norm(d_model)
        self.norm_3 = Norm(d_model)
        
        self.dropout_1 = nn.Dropout(dropout)
        self.dropout_2 = nn.Dropout(dropout)
        self.dropout_3 = nn.Dropout(dropout)
        
        self.attn_1 = MultiHeadAttention(heads, d_model, dropout=dropout)
        self.attn_2 = MultiHeadAttention(heads, d_model, dropout=dropout)
        self.ff = FeedForward(d_model, dropout=dropout)

    def forward(self, x, e_outputs, src_mask, trg_mask, Writedata):
        x2 = self.norm_1(x)
        x = x + self.dropout_1(self.attn_1(x2, x2, x2, trg_mask, Writedata))
        x2 = self.norm_2(x)
        x = x + self.dropout_2(self.attn_2(x2, e_outputs, e_outputs, \
        src_mask, Writedata))
        x2 = self.norm_3(x)
        x = x + self.dropout_3(self.ff(x2))
        return x