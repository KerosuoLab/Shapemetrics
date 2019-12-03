%% 1) Read in the Ilasti prediction maps for membrane (MEMB) and nuclei (NUCL)

% 1.1) PUT THE FILENAME OF ILASTIK PREDICTION MAP FOR MEMBRANE HERE:
%===================================================================%
ilastik_filename = 'imagename_Probabilities.h5';
%===================================================================%


% 1.3) USE h5read FUNCTION TO READ IN THE PROBABILITY MAPS:
ilastik_file = h5read(ilastik_filename,'/exported_data/'); % exported data is the "folder" in Ilastik project where the data is stored. Don't change this.
pred = squeeze(ilastik_file(2,:,:,:));                     % Originally 4D file, we want 3D
pred = permute(pred,[2,1,3]);                              % prevent the Ilastik axis swap by fixing them here
                                                                    
%% 2) VISUALIZE THE PREDICTION MAPS :
figure                                                              
imshow(pred)                                                        
title('Ilastik prediction map, 2D membrane ch')                     

%% 3) Read the original image in you used in Ilastik (the raw image)

% 3.1) PUT HERE THE ORIGINAL, RAW IMAGE FILENAME, MEMBRANE AND NUCLEUS:
%===================================================================%
imagename_2D = 'imagename.tif'; % membrane file name
%===================================================================%

% 3.2) READ IN THE IMAGES IN (FOR 2D IMAGES):
original_img_2D = imread(imagename_2D);                   

% 3.3) REUIRE THE PIXEL SIZE (SEE FROM FIJI)
pixelSize = 0.3352; % put your value here

%% 3.3) VISUALIZE THE JUST READ-IN IMAGES AS     

figure                                                              
imshow(original_img_2D)                                             
title('2D raw membrane image')                                      
                                                                    
%% 4) Take Ilastik prediction maps and threshold it: test with 5 diff. vals

%===================================================================%
% colors versus values:                                             %
%-------------------------------------------------------------------%
% o is black in output images                                       %
% 1 is white in output images                                       %
%-------------------------------------------------------------------%
% NOTE: the black and white are inverted between matlab and fiji!   %
% this means if you look output images in fiji, the colors are      %
% vice versa                                                        %
%===================================================================%

% 4.1) SET DIFFERENT THRESHOLD VALUE VARIABLES: (FIVE IN THIS CASE)

seg1 = pred>0.5; 
seg2 = pred>0.6;
seg3 = pred>0.7;
seg4 = pred>0.8;
seg5 = pred>0.9;

%% 4.2) VISUALIZE ALL THE THRESHOLD VALUES: MEMBRANE AND NUCLEUS
%===================================================================%
% the "figure" command in new lines creates own figure for each     %
% seg value instead of writing over the previous one                %
%===================================================================%
%===================================================================%
% threshold value >0.5 visualization:                               %
figure                                                              %
imshow(seg1)                                                   
title('seg1 MEMB (membrane ch), 2D, ilastik pred. map th >0.5')     
%% threshold >0.6 pixel values visualization:                       
figure                                                              
imshow(seg2)                                                   
title('seg2 MEMB (membrane ch), 2D, ilastik prediction map th >0.6')
%% threshold >0.7 pixel values visualization:                       
figure                                                              
imshow(seg3)                                                   
title('seg3 MEMB (membrane ch), 2D ilastik prediction map th >0.7') 
%% threshold >0.8 pixel values visualization:                       
figure                                                              
imshow(seg4)                                                   
title('seg4 MEMB (membrane ch), 2D, ilastik prediction map th >0.8')
%% threshold >0.9 pixel values visualization:                       
figure                                                              
imshow(seg5)                                                   
title('seg5 MEMB (membrane ch), 2D ilastik prediction map th>0.9')  
%% 5) Require now the minimum and maximum cell and nucleus size

% 5.1) SIZE THRESHOLD AND VISUALIZE ALL FIVE, MEMB and NUCL     
% seg1:                                                             
seg1 = bwareaopen(seg1,50);                 % min size      
seg1 = seg1 - bwareaopen(seg1,9000);        % max size         
figure                                                              
imshow(seg1)                                                   
title('2D, size th prediction map seg0 MEMB')                      

%% seg2:                                                            
seg2 = bwareaopen(seg2,50);                               
seg2 = seg2-bwareaopen(seg2,6000);  

figure                                                              
imshow(seg2)                                                   
title('2D, size th prediction map seg2 MEMB')                       
%% seg3:                                                            
seg3 = bwareaopen(seg3,50);                               
seg3 = seg3-bwareaopen(seg3,20000);                   
figure                                                              
imshow(seg3)                                                   
title('2D, size th prediction map seg3 MEMB')                       

%% seg4:                                                            
seg4 = bwareaopen(seg4,50);                               
seg4 = seg4-bwareaopen(seg4,2000);                   
figure                                                              
imshow(seg4)                                                   
title('2D, size th prediction map seg4 MEMB')                       
%% seg5:                                                            
seg5 = bwareaopen(seg5,50);                               
seg5 = seg5-bwareaopen(seg5,2000);                   
figure                                               
imshow(seg5)                                         
title('2D, size th prediction map seg5 MEMB')        

%===================================================================%
% NOTE: At this point, you choose the best pred threshold value. It %
% is by default the seg5 which gives the best results. Even higher  %
% thresholds can be tried if one wants.                             %
%===================================================================%

%% 7) CREATE THE LABEL MATRIX AND SIZE THRESHOLD IT:

seed   = imimposemin(original_img_2D,seg1);
Label  = watershed(seed);
Label2 = bwareaopen(Label,50);
Label2 = Label2 - bwareaopen(Label2,2000);
Final_Label = bwlabeln(Label2);   

% Add nice colormap for visualization:
Final_Label_r = label2rgb(Final_Label,'jet',[0.5,0.5,0.5]);

%% Visualize the segmentation

figure                % colorless
imshow(Label2)
title('Final Label')

figure                % with colors
imshow(Final_Label_r)
title('Final label, colored');
%% 7.5) EXTRACT THE SPATIAL INFORMATION OF EACH CELL:

stats = regionprops(Label_an,'All');
Areas = cat(1, stats.Area);
Number_Of_Cells = size(Areas,1);

%% 7.6) HISTOGRAM PLOTS FOR PARAMETER DISTRIBUTIONS
figure                                                             
hist(Areas*pixelValue*pixelValue,40) % the Areas matrix is in pixels squared, so multiplu with pixelValue^2                                         
title('Cell Areas, watershed (number of cells = )') % write here the number of cells
ylabel('Number of Cells with certain area')
xlabel('Area in microns squared (um)^2')
