## How to use
### Step 1: Change mode in Simulator/Makefile
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
* Notice that if you want to change the type of softmax simulation, Please modify the `bool STAR` in soft.cpp.
```c=
	bool STAR = false;	// <=== Modify When Change Structure ===>
```
### Step 2: Compile the NeuroSim & Check the shell file if needed
* After the setting, type `make` to compile the NeuroSim. if there is something wrong, type `make clean` and try again.
* Then check  data/Shell and find the corresponding .sh file.
    * TransSeg_Arch.sh $\Leftrightarrow$ mode=1
    * STAR_Arch.sh $\Leftrightarrow$ mode=2
    * STAR_Soft.sh & TransSeg_Soft.sh $\Leftrightarrow$ mode=3
* Mind that if you need to Simulate different numbers of Segment results, modify the TransSeg_Soft.sh's `argv[5] & argv[8]` to related number.

### Step 3: Run the Simulator
* After all, type the code as shown below, and the result will be displayed on your terminal.
```python=
python3 Simulate.py --mode={shell_name}
# Shell name list
# -- STAR_Arch
# -- STAR_Soft
# -- TransSeg_Arch
# -- TransSeg_Arch
```
