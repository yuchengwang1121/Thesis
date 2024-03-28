import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import numpy as np


class Norm(nn.Module):
    def __init__(self, d_model, eps = 1e-6):
        super().__init__()
    
        self.size = d_model
        
        # create two learnable parameters to calibrate normalisation
        self.alpha = nn.Parameter(torch.ones(self.size))
        self.bias = nn.Parameter(torch.zeros(self.size))
        
        self.eps = eps
    
    def forward(self, x):
        norm = self.alpha * (x - x.mean(dim=-1, keepdim=True)) \
        / (x.std(dim=-1, keepdim=True) + self.eps) + self.bias
        return norm

def attention(q, k, v, d_k, mask=None, dropout=None,Layer=None,  Encoder=None, Segnum=0):
    scores = torch.matmul(q, k.transpose(-2, -1)) /  math.sqrt(d_k)

    if mask is not None:
        mask = mask.unsqueeze(1)
        scores = scores.masked_fill(mask == 0, -1e9)
    if(Encoder and Segnum!=0):
        segment = scores.squeeze().detach().numpy()
        seg_size = int(segment.shape[0]/Segnum)

        merge_S= []
        SegSoft = np.empty((Segnum, Segnum, seg_size, seg_size))
        segment_R = np.split(segment, Segnum, axis=0)

        for i in range(len(segment_R)):
            for j in range(len(segment_R)):
                SegSoft[i][j] = F.softmax(torch.tensor(np.split(segment_R[i], Segnum, axis=1)[j]),dim=-1)
            merge_S.append(np.concatenate(SegSoft[i], axis=1))
        Res = np.concatenate(merge_S, axis=0)
        # np.savetxt('./data/Soft/soft_L'+str(Layer)+ '_'+str(Segnum) + '.csv', Res , delimiter=",",fmt='%10.5f')
        scores = torch.tensor(Res).to(v.dtype)
    else:
        scores = F.softmax(scores, dim=-1)
        # if(Encoder):
            # np.savetxt('./data/Soft/soft_L'+str(Layer)+ '_'+str(Segnum) + '.csv', scores.squeeze().detach().numpy() , delimiter=",",fmt='%10.5f')
        
    
    if dropout is not None:
        scores = dropout(scores)
    output = torch.matmul(scores, v)
    return output

class MultiHeadAttention(nn.Module):
    def __init__(self, heads, d_model, dropout = 0.1):
        super().__init__()
        
        self.d_model = d_model
        self.d_k = d_model // heads
        self.h = heads
        
        self.q_linear = nn.Linear(d_model, d_model)
        self.v_linear = nn.Linear(d_model, d_model)
        self.k_linear = nn.Linear(d_model, d_model)
        
        self.dropout = nn.Dropout(dropout)
        self.out = nn.Linear(d_model, d_model)
    
    def forward(self, q, k, v, mask=None, Writedata=None, Layer=0, Encoder=None, Segnum=0):
        w_filename = "./Weight/weight" + str(Layer)
        bs = q.size(0)
        # perform linear operation and split into N heads
        k = self.k_linear(k).view(bs, -1, self.h, self.d_k)
        q = self.q_linear(q).view(bs, -1, self.h, self.d_k)
        v = self.v_linear(v).view(bs, -1, self.h, self.d_k)
        if(Writedata):
            print("===> Writing the weight Q,K,V in ", w_filename)
            np.savetxt(w_filename + "_Q.csv", self.q_linear.weight.detach().numpy() , delimiter=",",fmt='%10.5f')
            np.savetxt(w_filename + "_K.csv", self.k_linear.weight.detach().numpy() , delimiter=",",fmt='%10.5f')
            np.savetxt(w_filename + "_V.csv", self.v_linear.weight.detach().numpy() , delimiter=",",fmt='%10.5f')
        # transpose to get dimensions bs * N * sl * d_model
        k = k.transpose(1,2)
        q = q.transpose(1,2)
        v = v.transpose(1,2)
        # print("===> The q.transpose is ", q.size())

        # calculate attention using function we will define next
        scores = attention(q, k, v, self.d_k, mask, self.dropout, Layer=Layer ,Encoder=Encoder, Segnum=Segnum)
        # concatenate heads and put through final linear layer
        concat = scores.transpose(1,2).contiguous()\
        .view(bs, -1, self.d_model)
        output = self.out(concat)
    
        return output

class FeedForward(nn.Module):
    def __init__(self, d_model, d_ff=2048, dropout = 0.1):
        super().__init__() 
    
        # We set d_ff as a default to 2048
        self.linear_1 = nn.Linear(d_model, d_ff)
        self.dropout = nn.Dropout(dropout)
        self.linear_2 = nn.Linear(d_ff, d_model)
    
    def forward(self, x):
        x = self.dropout(F.relu(self.linear_1(x)))
        x = self.linear_2(x)
        return x
