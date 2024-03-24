## How to use
### Step 1 : set the parameter
* The parameter is initialized as below, more detailed description can be found by using `-h`.
  * --Path: "../Data/"
  * --STAR: False
  * --TransSeg: False
  * --Seg: 4
  * --CalculMSE: False
### Step 2 : Generate the result of softmax
* To run the STAR softmax result
  ```python=
  python3 Softmax.py --STAR=true
  ```
* To run the TransSeg softmax result
  ```python=
  python3 Softmax.py --TransSeg=true
  ```
* To run the TransSeg softmax result with different number of segment
  ```python=
  python3 Softmax.py --TransSeg=True --Seg={number you want}
  ```
* To run the MSE
  ```python=
  python3 Softmax.py --CalculMSE=true
  ```
### Step 3 : Calculate the MSE between to file
  ```python=
  python3 Softmax.py --CalculMSE=true --src1={path of STAR_softmax} --src2=={path of TransSeg_softmax}
  ```
