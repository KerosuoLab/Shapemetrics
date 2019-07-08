# 3D-segmentation script
First, make sure that you are in the right folder (where all your files are so that matlab can download the files in without errors.
The folder can be changed from the small arrow on the top line of this window, where your current path is  showing. 
Run each section with command shift+enter. The most common error when trying to run the first section is being in the wrong folder,
which matlab states as follows: "the file xxx does not exist".
Only change names etc. that are between these lines.       
Anything outside the lines is not needed to change, and possible errors do not come from parts outside the lines.     

**Loading the Ilastik prediction map in**
```
ilastik_filename_MEMB = 'filename.h5';
ilastik_file_MEMB = h5read(ilastik_filename_MEMB,'/exported_data/');
pred_MEMB = squeeze(ilastik_file_MEMB(2,:,:,:));
pred_MEMB = permute(pred_MEMB,[2,1,3]);
```
Let's visualize prediction map as z-projection
```
figure                                                              
imshow(sum(pred_MEMB,3),[])                                         
title('Ilastik prediction map, z-projection, membrane ch')
```
![](images/tbud_ilastik_predmap_zproj.png)

