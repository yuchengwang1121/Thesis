## How to use
### Softmax.py
#### Step 1 : set the parameter
* The parameter is initialized as below, more detailed description can be found by using `-h`.
  * --WeightPath: "./data/Weight/"
  * --InputPath: "./data/Input/"
  * --ResultPath: "./data/Result/"
  * --Layer: 0
  * --STAR: False
  * --TransSeg: False
  * --Seg: 4
  * --CalculMSE: False
  * --src1: None
  * --src2: None
  * --GenCAMLUT: False
  * --CAM_size: 64
  * --write_CAM: False
  * --LUT_size: 16
  * --write_LUT: False
  * --Bin_Pre: 32
#### Step 2 : Run the needed function
* To run the STAR softmax result with different layer
  ```python=
  python3 Softmax.py --STAR=true --Layer={The Layer you want}
  ```
* To run the TransSeg softmax result
  ```python=
  python3 Softmax.py --TransSeg=true --Layer={The Layer you want}
  ```
* To run the TransSeg softmax result with different number of segment
  ```python=
  python3 Softmax.py --TransSeg=True --Seg={number you want} --Layer={The Layer you want}
  ```
* To run the MSE
  ```python=
  python3 Softmax.py --CalculMSE=true --src1={path of STAR_softmax} --src2=={path of TransSeg_softmax}
  ```
* To run the CAMLUT
  ```python=
  python3 Softmax.py --GenCAMLUT --write_CAM=true --write_LUT=true
  ```
* To run the CAMLUT with different size
  ```python=
  python3 Softmax.py --GenCAMLUT --CAM_size={your CAM size} --write_CAM=true --LUT_size={your LUT size} --write_LUT=true
  ```
* To run the CAMLUT with different precision
  ```python=
  python3 Softmax.py --GenCAMLUT --write_CAM=true --write_LUT=true --Bin_Pre={Precision you want}
  ```
### train.py
#### Step 1 : set the parameter
* The parameter is initialized as below.
  * -no_cuda: true
  * -SGDR: true
  * -epochs: 60
  * -d_model: 128
  * -n_layers: 5
  * -heads: 1
  * -dropout: 0.1
  * -batchsize: 1500
  * -printevery: 100
  * -lr:0.0001
  * -create_valset: true
  * -max_strlen: 80
  * -checkpoint: 0
  * -Writedata: False
#### Step 2 : Run the needed function
* Run first time
```python=
 # train.py -src_lang {src language you want} -trg_lang {target language you want} -src_data {src language txt} -trg_data {target language txt}
 python3 train.py -src_lang en_core_web_sm -trg_lang fr_core_news_sm -src_data data/english.txt -trg_data data/french.txt
```
* Run with pretrain model
```python=
# train.py -load_weights={The folder you save torch.dict}
python3 train.py -load_weights=Pretrain_weight -src_lang en_core_web_sm -trg_lang fr_core_news_sm -src_data data/english.txt -trg_data data/french.txt
```

### translate
#### Step 1 : set the parameter
* The parameter is initialized as below.
  * -k: 3
  * -max_len: 80
  * -d_model: 128
  * -n_layers: 5
  * -heads: 1
  * -dropout: 0.1
  * -no_cuda: true
  * -sim_ratio: 0.8
  * -test_times: 100
  * -Writedata: False
  * -Segnum: 0
#### Step 2 : Run the needed function
* Run
```python=
python3 translate.py -load_weights Pretrain_weight -src_lang en_core_web_sm -trg_lang fr_core_news_sm -Segnum={number you want}
```