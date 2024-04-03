## How to use
### Step 1: Change mode in Simulator/Makefile and compile if you modify .cpp file
* You can see the code in Makefile as shown below, please change the `mode` to compile the corresponding .cpp file.
```Makefile=
mode = 1    //<=== Only need to change this one
            // 1 for TransSeg Architecture
            // 2 for STAR Architecture
            // 3 for Softmax simulation
ifeq ($(mode),1)
	main=TransSeg
else ifeq ($(mode),2)
    main=STAR
else ifeq ($(mode),3)
	main=Soft
endif
```
### Step 2: Compile the NeuroSim & Check the shell file if needed
* After the setting, type `make` to compile the NeuroSim. if there is something wrong, type `make clean` and try again.
* Then hell name list are shown as below
    * Arch.sh $\Leftrightarrow$ mode=1, 2
    * Softmax.sh $\Leftrightarrow$ mode=3

### Step 3: Run the Simulator
* After all,pass arguments according to the parameters you need.And the result will be displayed on your terminal.
```python=
### Please note that if you intend to use Method "TransSeg", remember to provide the value for `--segnum`.
### Mode list
### 1. STAR_Arch
### 2. STAR_Soft
### 3. TransSeg_Arch
### 4. TransSeg_Soft

python3 Simulate.py --mode={shell_name} --segnum={numofsegment}
```
