# cell_tracking_metrics
After we obtain the tracking results from our cell tracking program, we should have two parts of files.
## First part: .txt file
The first part is a .txt file that contains each cell's birth, death and mitosis. 
![](https://github.com/hrlblab/cell_tracking_metrics/blob/master/asset/res_track.png)
### .txt file data structure
```
${AOGM_EVA}
   |
   └——testing_dataset
           └——01_RES
                └———res_track.txt          
```
## Second part: .tif file
The second part is consisted of several 16-bit .tif files that contain each cell's information corresponding to .txt file.
![](https://github.com/hrlblab/cell_tracking_metrics/blob/master/asset/mask.png)
### .tif file data structure
```
${AOGM_EVA}
   |
   └——testing_dataset
           └——01_RES
                └———maks000.tif-mask.0**.tif        
```
## Data structure
Before the evaluation, we need to organize all files into the structure below
```
${AOGM_EVA}
   |
   └—--—testing_dataset
             └——01_RES
                  └———maks000-mask0**.tif, res_track.txt          
             └——01_GT
                  └————SEG
                        └———man_seg000-man_seg0..                     
                  └————TRA
                        └———man_track000-man_track0.., man_track.txt  
```
## Evaluation
The software including introduction and routine can be downloaded here:[AOGM](https://drive.google.com/drive/folders/11tJ3qc2_D_R9ovCxJ6T8TbPkM1iLGXe1?usp=sharing)
We can use Terminal to execute the software in cell trakcing challenge with command line like below(Linux):
```
cd AOGM_EVA
./DETMeasure AOGM_EVA/testing_dataset 01 3
./SEGMeasure AOGM_EVA/testing_dataset 01 3
./TRAMeasure AOGM_EVA/testing_dataset 01 3
```
This software can be executed in Mac, Win and Linux. 
