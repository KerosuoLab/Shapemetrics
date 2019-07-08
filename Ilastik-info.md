# Ilastik instructions

- [A. Before Ilastik](#a-before-ilastik---first-steps)
- [B. Ilastik machine learning](#b-ilastik-machine-learning)

## A. Before Ilastik - first steps

The image file type used throughout our analysis is .tiff. However, raw microscope data rarely comes in .tif format by default. In addition, the stacks are often stored in a folder of individual 2D stacks instead of one 3D z-stack image. Convert the microscope file to .tif z-stack using the following steps:

**1.** Open Fiji/ImageJ and drag the raw microscope image folder of individual stacks to the program to open it.

**2.** Fiji/ImageJ will ask if the content of this folder will be converted into z-stack image. Press “yes”

**3.** In case of multiple channels merged together in one z-stack, separate the membrane channel by a) removing other stacks from “Image” > “Stacks” > “Tools” > “Slice remover” and entering the first and last slice number of the unnecessary slices, b) Only dragging the membrane stacks from the raw image folder in step 1.

**4.** At this point you should have z-stack image of membrane staining channel. It is recommended to use greyscale images in the analysis, so the image type should be changed to uint8 instead of RGB: “Image” > “Type” > “uint8”

**5.** Last step in FIJI/ImageJ is to save this z-stack as .tif: “File” > “Save as” > .tif.

##  B. Ilastik machine learning

The goal of this part is to perform the machine learning based segmentation for the greyscale z-stack image saved as .tif in the part A and export the resulting prediction map. Ilastik works with other formats and with single stacks / 2D images but the Matlab part of our analysis is designed for 3D .tif files only. 

**1.** Open Ilastik and choose “Pixel Calssification”

**2.** Program will ask you to save the teaching file to disk. Recommended (but not mandatory) is to save it in the same folder with your input file i.e. the z-stack image generated in the part A.

**3.** Pixel classification workflow is now opened. Upload your z-stack .tiff image file from “Add New” > “Add Separate Image(s)” and select the .tiff z-stack image file.

**4.** From “Feature Selection” we chose all the features in our pipeline. 

**5.** After selecting the features, move on to “Training” and add two labels by pressing “+Add Label”. With Label1 (red) mark the membrane and with Label2 (green) mark the background and the cell interiors. The result will be more accurate when the drawing is done precisely. After drawing, press “Live Update” to start the machine learning algorithm and visualize the results.

**6.** After the live update is done, the image will have red and green overlay. This indicates the result of machine learning segmentation. If you think the result could be better, draw more labels to the image to improve the result.

**7.** The uncertainty value can be checked by pressing the eye-icon next to “Uncertainty”. This will mark the pixels with uncertain segmentation result in turquoise. Reduce the amount of uncertainty by adding more labels.

**8.** After segmentation is finished, export the segmentation i.e. pixel value probability map from the section “Prediction export” and choosing “Export all”. This will generate .h5 file in the same folder as you original z-stack image.
