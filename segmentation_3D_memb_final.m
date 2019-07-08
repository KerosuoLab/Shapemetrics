%% 1) Read in the Ilastik prediction maps for membrane (MEMB) and nuclei (NUCL)

%___________________________________________________________________%
                                                                    %
% First, make sure that you are in the right folder (where all your %
% files are so that matlab can download the files in without errors.%
% The folder can be changed from the small arrow on the top line of %
% this window, where your current path is  showing. Run each section%
% with command shift+enter.                                         %
% The most common error when trying to run the first section is     %
% being in the wrong folder, which matlab states as follows:        %
% "the file xxx does not exist".                                    %
                                                                    %
% Only change names etc. that are between these "====" lines.       %
% Anything outside the "====" lines is not needed to change, and    %
% possible errors do not come from parts outside the lines.         %
                                                                    %
% The text between "_____" lines are info that you should read      %
% before running that section.                                      %
%___________________________________________________________________%


% 1.1) PUT THE FILENAME OF ILASTIK PREDICTION MAP FOR MEMBRANE HERE:
%===================================================================%
ilastik_filename_MEMB = 'filename.h5'; % ilastik prediction map here
%===================================================================%

% 1.2) USE h5read FUNCTION TO READ IN THE ILASTIK PROBABILITY MAP:
ilastik_file_MEMB = h5read(ilastik_filename_MEMB,'/exported_data/'); % exported data is the "folder" in Ilastik project where the data is stored. Don't change this.
pred_MEMB = squeeze(ilastik_file_MEMB(2,:,:,:));                     % Originally 4D file, we want 3D
pred_MEMB = permute(pred_MEMB,[2,1,3]);                              % prevent the Ilastik axis swap by fixing them here

% 1.3) VISUALIZE THE SUMMED Z-PROJECTIONS OF PREDICTION MAPS :
figure                                                              
imshow(sum(pred_MEMB,3),[])                                         
title('Ilastik prediction map, z-projection, membrane ch')          
%% 2) Read the original image in you used in Ilastik (the raw image)

%___________________________________________________________________%
% The original image must be read in inside a loop since 3D images  %
% are just bach of 2D images in stack in z-direction. We loop thus  %
% through the z-direction.                                          %
                                                                    %
% NOTE: the figure visualizations are separated in different        % 
% sections to ensure that the titles are in correct images.         %
%___________________________________________________________________%

% 2.1) PUT HERE THE ORIGINAL, RAW IMAGE FILENAME OF MEMBRANE:
%===================================================================%
imagename_MEMB = 'imagename.tif';   % membrane file name
%===================================================================%

% 2.2) READ IN THE IMAGE INSIDE THE LOOP (LOOP FOR 3D IMAGES):
original_img_MEMB = 0*pred_MEMB;                                    
for z = 1 : size(pred_MEMB,3)                % We loop through the z-direction: size(pred,3) means the length of z-axis in pred-file                          
    temp = imread(imagename_MEMB,z);         % temp file is each 2D stack we read in and write over in each run of the loop
    original_img_MEMB(:,:,z) = temp(:,:,1);  % before writing over the 2D stack, we save it here in 3D matrix as one stack
end                                                                 

% 2.3) VISUALIZE THE JUST READ-IN IMAGE AS SUMMED MAX-INTENSITY and SUMMED ZPROJECTIONS:     
figure                                                              
imshow(max(original_img_MEMB,[],3),[])                             
title('max-projection, raw membrane image')                         
figure
imshow(sum(original_img_MEMB,3),[])
title('z-projection, raw membrane image')
%% 3) THE BLURRED VERSIONS OF RAW IMAGES TO CONNECT ANY GAPS
% 3.1) USE STREL 3D TO BLUR AND THUS CONNECT ANY GAPS IN THE MEMBRANE:
img_blur_MEMB = imdilate(original_img_MEMB,strel3D('sphere',3));    

% 3.2) VISUALIZE THE BLURRED VERSIONS AS Z-PROJECTIONS:
figure                                                              
imshow(sum(img_blur_MEMB,3),[])                                     
title('Blurred version of membrane image, using strel3D function')  
%% 4) Take Ilastik prediction map and threshold it: test with 5 diff. vals

%___________________________________________________________________%
% colors versus values:                                             %
%-------------------------------------------------------------------%
% o is black in output images                                       %
% 1 is white in output images                                       %
%-------------------------------------------------------------------%
% NOTE: the black and white are inverted between matlab and fiji!   %
% this means if you look output images in fiji, the colors are      %
% vice versa                                                        %
%___________________________________________________________________%

% 4.1) SET DIFFERENT THRESHOLD VALUE VARIABLES: (SIX IN THIS CASE)

seg1_MEMB = pred_MEMB>0.5; 
seg2_MEMB = pred_MEMB>0.6;
seg3_MEMB = pred_MEMB>0.7;
seg4_MEMB = pred_MEMB>0.8;
seg5_MEMB = pred_MEMB>0.9;
seg6_MEMB = pred_MEMB>0.95;

%% 4.2) VISUALIZE ALL THE THRESHOLD VALUES:
% threshold value >0.5 visualization:                               
figure                                                              
imshow(sum(single(seg1_MEMB),3),[]);                                
title('seg1 MEMB (membrane ch), ilastik pred. map th >0.5')         
%% threshold >0.6 pixel values visualization:                       
figure                                                              
imshow(sum(single(seg2_MEMB),3),[]);                                
title('seg2 MEMB (membrane ch), ilastik prediction map th >0.6')    
%% threshold >0.7 pixel values visualization:                       
figure                                                              
imshow(sum(single(seg3_MEMB),3),[]);                                
title('seg3 MEMB (membrane ch), ilastik prediction map th >0.7')    
%% threshold >0.8 pixel values visualization:                       
figure                                                              
imshow(sum(single(seg4_MEMB),3),[]);                                
title('seg4 MEMB (membrane ch), ilastik prediction map th >0.8')    
%% threshold >0.9 pixel values visualization:                       
figure                                                              
imshow(sum(single(seg5_MEMB),3),[]);                                
title('seg5 MEMB (membrane ch), ilastik prediction map th>0.9')     
%%                                                                  
figure                                                                                 
imshow(sum(single(seg6_MEMB),3),[]);                                
title('seg6 MEMB (membrane ch), ilastik prediction map th>0.95')    
%% 4.3) Require now the minimum and maximum cell size - visualize them
 
%___________________________________________________________________%
% this is the tricky part. The size threshold can be changed manually:
% The whole point is to rmove the background and leave the separated cells.

% run with the value that is already here first, after which if you see nothing,
% change the max% value to be smaller (remove one zero). 
% If you see the background (big 
% white block around the sample) make the max size threshold bigger by one
% zero. 

% NOTE: The seg1 very likely doesn't give anything because the probability
% map threshold is so low. In most cases, only the seg4-seg6 give
% you good images, so do not worry if the three first thresholds give you
% only balck or only background.

% NOTE: if you change the threshold to be smaller, run the part 4.1. 
% before running the changed threshold again, otherwise there is no change
% to the last attempt. After becoming familiar with this part, it is
% usually worth of running only se4-se6 parts.
%___________________________________________________________________%

% seg1:
%===================================================================%
seg1_MEMB = bwareaopen(seg1_MEMB,100);                 % min size   
seg1_MEMB = seg1_MEMB-bwareaopen(seg1_MEMB,200000);    % max size   
%===================================================================%

figure                                                             
imshow(sum(single(seg1_MEMB),3),[]);                               
title('z-projection of size th prediction map seg1 MEMB, pixel values > 0.5')          
%% seg2:                                     
%===================================================================%
seg2_MEMB = bwareaopen(seg2_MEMB,100);                                
seg2_MEMB = seg2_MEMB-bwareaopen(seg2_MEMB,200000);      
%===================================================================%
figure                                                              
imshow(sum(single(seg2_MEMB),3),[]);                                
title('z-projection of size th prediction map seg2 MEMB, pixel values > 0.6')           
%% seg3:                                             
%===================================================================%
seg3_MEMB = bwareaopen(seg3_MEMB,100);                               
seg3_MEMB = seg3_MEMB-bwareaopen(seg3_MEMB,200000);      
%===================================================================%
figure                                                              
imshow(sum(single(seg3_MEMB),3),[]);                                
title('z-projection of size th prediction map seg3 MEMB,, pixel values > 0.7')           
%% seg4:                                  
%===================================================================%
seg4_MEMB = bwareaopen(seg4_MEMB,100);                               
seg4_MEMB = seg4_MEMB-bwareaopen(seg4_MEMB,200000);      
%===================================================================%
figure                                                              
imshow(sum(single(seg4_MEMB),3),[]);                                
title('z-projection of size th prediction map seg4 MEMB, pixel values > 0.8')           
%% seg5:                                   
%===================================================================%
seg5_MEMB = bwareaopen(seg5_MEMB,100);                               
seg5_MEMB = seg5_MEMB-bwareaopen(seg5_MEMB,200000);   
%===================================================================%
figure                                                              
imshow(sum(single(seg5_MEMB),3),[]);                                
title('z-projection of size th prediction map seg5 MEMB,, pixel values > 0.9')           
%% seg6            
%===================================================================%
seg6_MEMB = bwareaopen(seg6_MEMB,900);                               
seg6_MEMB = seg6_MEMB - bwareaopen(seg6_MEMB,100000);   
%===================================================================%
figure                                                              
imshow(sum(single(seg6_MEMB),3),[]);                                
title('z-projection of the size th prediction map seg6 MEMB, pixel values >0.95')
%% 5) CHOOSE THE BEST THRESHOLD VALUE TO USE BASED ON HOW MANY CELLS ARE STILL THERE AFTER BACKGROUND IS GONE
%___________________________________________________________________%
% NOTE: At this point, you choose the best pred threshold value. It %
% is by default the seg6 which is noticed to give best results      %
% (but this might varie depending on the sample).                   %
% If you want to use another seg, write it below                    %
%___________________________________________________________________%

% PUT HERE THE "SEG" THAT YOU WANT TO USE: 
%===================================================================%
seg_final = seg6_MEMB; % alternaitves can be seg5_MEMB etc.
%===================================================================%
%% 6) Create the seed for waterhsed algorithm by using the blurred image and chosen "seg"

seed_MEMB = imimposemin(img_blur_MEMB,seg_final);                    
%% 7) Create the label matrix using watershed and extract spatial info from the label

% 7.1) CREATE THE LABEL MATRIX:
Label_MEMB = watershed(seed_MEMB);                               

% 7.2) SIZE THRESHOLD THE LABEL MATRIX TO EXTRACT THE BACKGROUND VOLUMES:
%___________________________________________________________________%
% This size thresholding can be manually changed with your opinion how good
% does the label matrix look: IT CAN BE VIEWED IN MATLAB APP "VOLUME
% VIEWER" as 3d rendition

% First keep this size threshold same as the seg threshold
% If the 3D rendition does not look good, you ca nchange the size threshold
% below
%___________________________________________________________________%

%===================================================================%
Label2_MEMB = bwareaopen(Label_MEMB,200);                    % min                   
Label2_MEMB = Label2_MEMB - bwareaopen(Label_MEMB,200000);   % max   
%===================================================================%

% 7.3) LABEL THE SIZE THRESHOLDED MATRIX TO GET INFO FROM INDIVIDUAL CELLS:
Final_Label_MEMB = bwlabeln(Label2_MEMB);                          

%% 8) SAVE THE LABELS AND SEGMENTATION VISUALIZATION TO DISK
% 8.1) save the final label as tiff z-stack to disk, look at it in fiji:
for z = 1 : size(Final_Label_MEMB,3)                                
    imwrite(Final_Label_MEMB(:,:,z),'Final_Label_membrane.tif','compression','none','writemode','append');
end                                                                
% 8.2) Save segmentation borders as tiff to disk, look at them in fiji:  
for z = 1 : size(original_img_MEMB,3)                               
    temp        = zeros(size(original_img_MEMB,1),size(original_img_MEMB,2),3,'uint8');
    per         = Final_Label_MEMB(:,:,z) == 0;                             
    temp(:,:,1) = original_img_MEMB(:,:,z);      % original image on the back on red
    temp(:,:,3) = uint8(per)*100;                % Label borders on blue on top   
    imwrite(temp,'Segmentation_borders_membrane.tif','compression','none','WriteMode','append'); 
end                                                                 
                                                           
%% 9) EXTRACT THE SPATIAL INFORMATION OF EACH CELL:

stats_MEMB = regionprops3(Final_Label_MEMB,'all');

CellVolumes          = stats_MEMB.Volume;                   % cell volumes
CellSurfaceAreas     = stats_MEMB.SurfaceArea;              % cell surface area
CellCentroids        = stats_MEMB.Centroid;                 % cell centroids
CellVolSurfAreaRatio = CellVolumes./CellSurfaceAreas;       % cell volume-surface area ratio
CellEllipticity      = (stats_MEMB.PrincipalAxisLength(:,1) - stats_MEMB.PrincipalAxisLength(:,3))./(stats_MEMB.PrincipalAxisLength(:,1)); % ellipticity
LongestAxis          = stats_MEMB.PrincipalAxisLength(:,1); % length of the longest axis (diameter of longest axis)
CellElongation       = LongestAxis./((stats_MEMB.PrincipalAxisLength(:,2).*stats_MEMB.PrincipalAxisLength(:,3))./2); % cell elongation: longets axis divided by the average of thwo shortest
NumberOfCells        = size(stats_MEMB.Volume,1);           % number of cells

%% 9.2) HISTOGRAM PLOT FOR VOLUME DISTRIBUTION
figure                                                              
hist(CellVolumes,100)

title('Cell volumes, number of cells = ')                       
ylabel('Number of Cells with certain volume')
xlabel('Cell Volume in voxels')
%% 10) Preparation for visualization of spatial parameter statistics

% 10.1) All the parameter values
stats_matrix_MEMB_all      = zeros(NumberOfCells,5);                    
stats_matrix_MEMB_all(:,1) = CellVolumes;                % parameter 1: volume              
stats_matrix_MEMB_all(:,2) = CellVolSurfAreaRatio;       % parameter 2: vol-surfArea ratio              
stats_matrix_MEMB_all(:,3) = CellEllipticity;            % parameter 3: ellipticity          
stats_matrix_MEMB_all(:,4) = CellElongation;             % parameter 4: elongation
stats_matrix_MEMB_all(:,5) = LongestAxis;                % parameter 5: the length of the longest axis
zscored_MEMB               = zscore(stats_matrix_MEMB_all);  

% 10.2) Knock out some parameters, six possibilities we may look:
stats_matrix_MEMB_1thru2 = stats_matrix_MEMB_all(:,1:2); % parameters 1 and 2
stats_matrix_MEMB_1thru3 = stats_matrix_MEMB_all(:,1:3); % parameters 1, 2 and 3   
stats_matrix_MEMB_1thru4 = stats_matrix_MEMB_all(:,1:4); % parameters 1, 2, 3 and 4
stats_matrix_MEMB_2thru5 = stats_matrix_MEMB_all(:,2:5); % parameters 2, 3, 4 and 5
stats_matrix_MEMB_3thru5 = stats_matrix_MEMB_all(:,3:5); % parameters 3, 4 and 5
stats_matrix_MEMB_4thru5 = stats_matrix_MEMB_all(:,4:5); % parameters 4 and 5
stats_matrix_MEMB_1and3and4 = stats_matrix_MEMB_all(:,[1,3,4]); % parameters 1, 3 and 4
stats_matrix_MEMB_1and4 = stats_matrix_MEMB_all(:,[1,4]); % parameters 1 and 4
stats_matrix_MEMB_1and3 = stats_matrix_MEMB_all(:,[1,3]); % parameters 1 and 3

% 10.3) zscore all of these partial number of spatial parameters:
zscored_MEMB_1thru2 = zscore(stats_matrix_MEMB_1thru2); 
zscored_MEMB_1thru3 = zscore(stats_matrix_MEMB_1thru3);
zscored_MEMB_1thru4 = zscore(stats_matrix_MEMB_1thru4);
zscored_MEMB_2thru5 = zscore(stats_matrix_MEMB_2thru5);
zscored_MEMB_3thru5 = zscore(stats_matrix_MEMB_3thru5);
zscored_MEMB_4thru5 = zscore(stats_matrix_MEMB_4thru5);
zscored_MEMB_1and3and4 = zscore(stats_matrix_MEMB_1and3and4);
zscored_MEMB_1and4 = zscore(stats_matrix_MEMB_1and4);
zscored_MEMB_1and3 = zscore(stats_matrix_MEMB_1and3);

%% 11) Save the stats-matrices into own folders for each sample: save here the label matrix and stats-matrix for the visualization purposes
% In addition, saving these will ensure you to be able to keep the data
% stored, not overwrited when analysing the next sample. way of doing this
% is: save('name to save as','name in the workspace')
% NOTE: we save them with same name as in the workspace:

save('Final_Label_MEMB','Final_Label_MEMB')
save('stats_MEMB','stats_MEMB')
save('stats_matrix_MEMB_all','stats_matrix_MEMB_all')

%% 12) VISUALIZATION: Clustergaram heatmaps of each single sample image: this is done once for each sample type 

% 12.1) name the parameters for heatmap in the order we have set above:
parameters_MEMB      = {'Cell Volume','Cell Volume/Surface ratio','Cell Ellipticity','Cell Elongation','Longest Axis'};
parameters_1thru2    = {'Cell Volume','Cell Volume/Surface ratio'};
parameters_1thru3    = {'Cell Volume','Cell Volume/Surface ratio','Cell Ellipticity'};
parameters_1thru4    = {'Cell Volume','Cell Volume/Surface ratio','Cell Ellipticity','Cell Elongation'};
parameters_2thru5    = {'Cell Volume/Surface ratio','Cell Ellipticity','Cell Elongation','Longest Axis'};
parameters_3thru5    = {'Cell Ellipticity','Cell Elongation','Longest Axis'};
parameters_4thru5    = {'Cell Elongation','Longest Axis'};
parameters_1and3and4 = {'Cell Volume','Cell Ellipticity','Cell Elongation'};
parameters_1and4     = {'Cell Volume','Cell Elongation'};
parameters_1and3     = {'Cell Volume','Cell Ellipticity'};

%%
% 12.2) create heatmaps from all of these parameter-mixes
heatm_MEMB_all       = clustergram(zscored_MEMB','RowLabels',parameters_MEMB','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1thru2    = clustergram(zscored_MEMB_1thru2','RowLabels',parameters_1thru2','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1thru3    = clustergram(zscored_MEMB_1thru3','RowLabels',parameters_1thru3','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1thru4    = clustergram(zscored_MEMB_1thru4','RowLabels',parameters_1thru4','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_2thru5    = clustergram(zscored_MEMB_2thru5','RowLabels',parameters_2thru5','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_3thru5    = clustergram(zscored_MEMB_3thru5','RowLabels',parameters_3thru5','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_4thru5    = clustergram(zscored_MEMB_4thru5','RowLabels',parameters_4thru5','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1and3and4 = clustergram(zscored_MEMB_1and3and4','RowLabels',parameters_1and3and4','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1and4     = clustergram(zscored_MEMB_1and4','RowLabels',parameters_1and4','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1and3     = clustergram(zscored_MEMB_1and3','RowLabels',parameters_1and3','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);

%% 12.3) Visualization info: TELL WHAT HEATMAP DO YOU WANT TO VISUALIZE

% HERE WRITE THE CORRECT HEATMAP INFO YOU WANT TO VISUALIZE:
%===================================================================%
% 1) what heatmap? (heatm = one out of the 7 heatmaps you have)add "all" or
% 1thru4 or 2thru4 etc.. to the end
heatm_to_visualize = heatm_MEMB_1and3;

% 2) what is the corresponding stat_matrix of this heatmap? add "all" or
% 1thru4 or 2thru4 etc.. to the end
stats_matrix_to_visualize = stats_matrix_MEMB_1and3;

% 3) what are the brances in this heatmap that you want to visualise?
% Choose from your heatmap!
branches = [3676 3661 3679 3656]; 

% 4) what is going to be the name of the colored .tif image of sub labels?
name = 'Stacks_heatm10_all.tif';
%===================================================================%
%% 12.4) Visalization itself: 
% 12.5) create an empty struct you are then filling in with parameter info needed
% to visualise the heatmap groups in original image:

xp_MEMB = struct('stats_MEMB_all',[],'SpatParamVals_MEMB',[],'CellIdentities',[],'Centroid',[],...
    'SpatParamVals_MEMB_len',[]);  

% 12.6) Fill this struct with info you told above
xp_MEMB.stats_MEMB_all         = load('stats_MEMB');
xp_MEMB.Centroid               = xp_MEMB.stats_MEMB_all.stats_MEMB.Centroid;
xp_MEMB.SpatParamVals_MEMB     = stats_matrix_to_visualize;
xp_MEMB.CellIdentities         = find(xp_MEMB.stats_MEMB_all.stats_MEMB.Volume);
xp_MEMB.SpatParamVals_t        = xp_MEMB.SpatParamVals_MEMB';
xp_MEMB.SpatParamVals_MEMB_len = size(xp_MEMB.SpatParamVals_MEMB,1);


% 12.7) Get the heatmap you stated above and assign colormap for the
% picture: note that these are not only colours you can use.

get(heatm_to_visualize)
cmp = jet(length(branches));
cmp(1,:)  = [0.1,0.9,1];   % light blue    
cmp(2,:)  = [0.8,0.6,0.9]; % lilac
cmp(3,:)  = [0,0.8,0.8];   % cyan   
cmp(4,:)  = [0.8,0.2,0.4]; % red pink   
cmp(5,:)  = [0.9,1,0];     % yellow
cmp(6,:)  = [0.1,0.4,1];   % dark blue     
cmp(7,:)  = [0.6,0.9,0.4]; % green
cmp(8,:)  = [0.9,0,1];     % purple 
cmp(9,:)  = [1,0.6,0];     % orange   
cmp(10,:) = [0.9,0.3,0.7]; % pink                   
cmp(11,:) = [1,0.7,0.2];   % diff. orange     
cmp(12,:) = [0.9,0.1,0.4]; % dif dark red     
cmp(13,:) = [0.5,1,0.6];   % different bright green    
cmp(14,:) = [0.9,0.4,0.5]; % baby pink       
cmp(15,:) = [0.8,0.2,0.4]; % red pink
cmp(16,:) = [0.1,0.9,1];   % light blue  

% The counter is like an arrow pointing that from which colormap color do
% you want to start. If counter = 1, means your first colour to use is
% cmp(1,:) which is currently light blue. 
%===================================================================%
counter=1;
%===================================================================%


sub_cluster_MEMB = struct('cells_of_interest',[]);
CellIdentities_MEMB = cat(1,xp_MEMB.CellIdentities);
range_MEMB = [0;cumsum(cat(1,xp_MEMB.SpatParamVals_MEMB_len))];
figure
imshow(max(original_img_MEMB,[],3),[]) 
hold on
for n = branches 
        group_of_interest_MEMB = clusterGroup(heatm_to_visualize, n, 'col');
        Col_Labels_MEMB = group_of_interest_MEMB.ColumnLabels; % here we have extracted some column labels from the clustergram
        Double_Labels_MEMB = cell(0);
        for i = 1 : length(Col_Labels_MEMB)
            Double_Labels_MEMB{i} = str2double(Col_Labels_MEMB{i});
        end
        goi_MEMB = cell2mat(Double_Labels_MEMB);
        cells_of_Interest_MEMB = CellIdentities_MEMB(goi_MEMB);
        sub_cluster_MEMB(counter).cells_of_interest = cells_of_Interest_MEMB;
    plot(xp_MEMB.stats_MEMB_all.stats_MEMB.Centroid(sub_cluster_MEMB(counter).cells_of_interest,1),...
         xp_MEMB.stats_MEMB_all.stats_MEMB.Centroid(sub_cluster_MEMB(counter).cells_of_interest,2),'*','color',cmp(counter,:),'LineWidth',3); 
    counter = counter+1;
end
%hold off
sub_cluster_pruned_MEMB = sub_cluster_MEMB;
for i = 1 : length(sub_cluster_MEMB)
    template = sub_cluster_MEMB(i).cells_of_interest;
    for j = 1 : length(sub_cluster_MEMB)
        pattern = sub_cluster_MEMB(j).cells_of_interest;
        if isempty(setdiff(pattern,template)) &&(i~=j)
           sub_cluster_pruned_MEMB(i).cells_of_interest = setdiff(template,pattern); 
        end
    end
    cells_of_Interest_MEMB = sub_cluster_pruned_MEMB(i).cells_of_interest; 
    in_MEMB = intersect(find(cells_of_Interest_MEMB>range_MEMB(1)),find(cells_of_Interest_MEMB<range_MEMB(2)));      
    template = xp_MEMB.CellIdentities; 
    pattern = cells_of_Interest_MEMB(in_MEMB);
    D = pdist2(template,pattern);
    in_xp_num = find(min(D,[],2)==0);
    id_MEMB = xp_MEMB.CellIdentities(in_xp_num);
    sub_cluster_pruned_MEMB(i).cells_of_interest = id_MEMB; 
end
Label_sub_MEMB = 0*Final_Label_MEMB;   
for c = 1 : length(sub_cluster_pruned_MEMB)
    cells_of_Interest_MEMB = sub_cluster_pruned_MEMB(c).cells_of_interest;
    for i = 1 : length(cells_of_Interest_MEMB) 
        Label_sub_MEMB(Final_Label_MEMB == cells_of_Interest_MEMB(i)) = c;
    end
end

for z = 1 : size(Label_sub_MEMB,3)
    temp  = zeros(size(Label_sub_MEMB,1),size(Label_sub_MEMB,2),3,'uint8');
    for c = 1 : length(sub_cluster_pruned_MEMB)
        tmp = Label_sub_MEMB(:,:,z);
        tmp(tmp~=c) = 0;
        tmp = tmp>0;
        temp(:,:,1) = uint8(cmp(c,1).*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(cmp(c,2).*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(cmp(c,3).*255*double(tmp))+temp(:,:,3);
    end
    imwrite(temp,name,'tiff','Compression','none','WriteMode','append');
end

%% STREL3D script (added from additional script strel3D.m) !!! NO NEED TO RUN THIS !!!
% this only a function, can be ignored by the user. Credit of this is for
% this function is for the SGA code makers

function se = strel3D(shape, size)
 N = size;   
    if strcmp(shape, 'sphere')        
        se = false([2*N+1 2*N+1 2*N+1]);
        [X,Y,Z] = meshgrid(-N:N, -N:N, -N:N);
        se(X.^2 + Y.^2 + Z.^2 <= N^2) = 1;
    elseif strcmp(shape, 'cube')
        se = true([2*N+1 2*N+1 2*N+1]);
    else 
        error('strel type not recognized');
    end
end
