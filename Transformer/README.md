## How to use
### Step 1 : set the parameter
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
### Step 2 : Run the needed function
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
