%% 1) Read in the Ilastik prediction maps for membrane 

% 1.1) PUT THE FILENAME OF ILASTIK PREDICTION MAP FOR MEMBRANE HERE:
%===================================================================%
ilastik_filename = 'nucl_trial_bin2_Probabilities.h5'; % ilastik prediction map here
%===================================================================%

% 1.2) USE h5read FUNCTION TO READ IN THE ILASTIK PROBABILITY MAP:
ilastik_file = h5read(ilastik_filename,'/exported_data/'); % exported data is the "folder" in Ilastik project where the data is stored. Don't change this.
pred = squeeze(ilastik_file(2,:,:,:));                     % Originally 4D file, we want 3D
pred = permute(pred,[2,1,3]);                              % prevent the Ilastik axis swap by fixing them here

% 1.3) VISUALIZE THE SUMMED Z-PROJECTIONS OF PREDICTION MAPS :
figure                                                              
imshow(sum(pred,3),[])                                         
title('Ilastik prediction map, z-projection')          
%% 2) Read the original image in you used in Ilastik (the raw image)

% 2.1) PUT HERE THE ORIGINAL, RAW IMAGE FILENAME OF MEMBRANE:
%===================================================================%
imagename = 'nucl_trial_bin2.tif';   % membrane file name
%===================================================================%

% 2.2) READ IN THE IMAGE INSIDE THE LOOP (LOOP FOR 3D IMAGES):
original_img = 0*pred;                                    
for z = 1 : size(pred,3)                % We loop through the z-direction: size(pred,3) means the length of z-axis in pred-file                          
    temp = imread(imagename,z);         % temp file is each 2D stack we read in and write over in each run of the loop
    original_img(:,:,z) = temp(:,:,1);  % before writing over the 2D stack, we save it here in 3D matrix as one stack
end                                                                 

% 2.3) VISUALIZE THE JUST READ-IN IMAGE AS SUMMED MAX-INTENSITY and SUMMED ZPROJECTIONS:     
figure                                                              
imshow(max(original_img,[],3),[])                             
title('max-projection, raw image')                         
figure
imshow(sum(original_img,3),[])
title('z-projection, raw image')
%% 3) THE BLURRED VERSIONS OF RAW IMAGES TO CONNECT ANY GAPS

% 3.1) USE STREL 3D TO BLUR AND THUS CONNECT ANY GAPS IN THE MEMBRANE:
img_blur = imdilate(original_img,strel3D('sphere',3));    

% 3.2) VISUALIZE THE BLURRED VERSIONS AS Z-PROJECTIONS:
figure                                                              
imshow(sum(img_blur,3),[])                                     
title('Blurred version of original image, using strel3D function')  
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

seg1 = pred>0.7;
seg2 = pred>0.8;
seg3 = pred>0.9;
seg4 = pred>0.95;

%% 4.2) VISUALIZE ALL THE THRESHOLD VALUES:
% threshold >0.7 pixel values visualization:                       
figure                                                              
imshow(sum(single(seg1),3),[]);                                
title('seg1, ilastik prediction map th >0.7')    
%% threshold >0.8 pixel values visualization:                       
figure                                                              
imshow(sum(single(seg2),3),[]);                                
title('seg2, ilastik prediction map th >0.8')    
%% threshold >0.9 pixel values visualization:                       
figure                                                              
imshow(sum(single(seg3),3),[]);                                
title('seg3, ilastik prediction map th>0.9')     
%%                                                                  
figure                                                                                 
imshow(sum(single(seg4),3),[]);                                
title('seg4, ilastik prediction map th>0.95')    
%% 4.3) Require now the minimum and maximum cell size - visualize them
          
% seg1:                                             
%===================================================================%
seg1 = bwareaopen(seg1,100);                               
seg1 = seg1-bwareaopen(seg1,15000);      
%===================================================================%
figure                                                              
imshow(sum(single(seg1),3),[]);                                
title('z-projection of size th prediction map seg1, pixel values > 0.7')           
%% seg4:                                  
%===================================================================%
seg2 = bwareaopen(seg2,100);                               
seg2 = seg2-bwareaopen(seg2,15000);      
%===================================================================%
figure                                                              
imshow(sum(single(seg2),3),[]);                                
title('z-projection of size th prediction map seg2, pixel values > 0.8')           
%% seg5:                                   
%===================================================================%
seg3 = bwareaopen(seg3,100);                               
seg3 = seg3-bwareaopen(seg3,15000);   
%===================================================================%
figure                                                              
imshow(sum(single(seg3),3),[]);                                
title('z-projection of size th prediction map seg3, pixel values > 0.9')           
%% seg6            
%===================================================================%
seg4 = bwareaopen(seg4,100);                               
seg4 = seg4 - bwareaopen(seg4,15000);   
%===================================================================%
figure                                                              
imshow(sum(single(seg4),3),[]);                                
title('z-projection of the size th prediction map seg4 MEMB, pixel values >0.95')
%% 5) CHOOSE THE BEST THRESHOLD VALUE TO USE BASED ON HOW MANY CELLS ARE STILL THERE AFTER BACKGROUND IS GONE

% PUT HERE THE "SEG" THAT YOU WANT TO USE: 
%===================================================================%
seg_final = seg2; 
%===================================================================%
%% 6) Create the seed for waterhsed algorithm by using the blurred image and chosen "seg"

seed = imimposemin(original_img,seg_final);                    
%% 7) Create the label matrix using watershed and extract spatial info from the label

% 7.1) CREATE THE LABEL MATRIX:
Label = watershed(seed);                               

% 7.2) SIZE THRESHOLD THE LABEL MATRIX TO EXTRACT THE BACKGROUND VOLUMES:

%===================================================================%
Label2 = bwareaopen(Label,600);             % min                   
Label2 = Label2 - bwareaopen(Label,15000);  % max   
%===================================================================%

% 7.3) LABEL THE SIZE THRESHOLDED MATRIX TO GET INFO FROM INDIVIDUAL CELLS:
Final_Label = bwlabeln(Label2);                          

%% 8) SAVE THE LABELS AND SEGMENTATION VISUALIZATION TO DISK

% 8.1) save the final label as tiff z-stack to disk, look at it in fiji:
for z = 1 : size(Final_Label,3)                                
    imwrite(Final_Label(:,:,z),'Final_Label_membrane.tif','compression','none','writemode','append');
end                                                                
% 8.2) Save segmentation borders as tiff to disk, look at them in fiji:  
for z = 1 : size(original_img,3)                               
    temp        = zeros(size(original_img,1),size(original_img,2),3,'uint8');
    per         = Final_Label(:,:,z) == 0;                             
    temp(:,:,1) = original_img(:,:,z);      % original image on the back on red
    temp(:,:,3) = uint8(per)*100;           % Label borders on blue on top   
    imwrite(temp,'Segmentation_borders_membrane.tif','compression','none','WriteMode','append'); 
end                                                                 
                                                           
%% 9) EXTRACT THE SPATIAL INFORMATION OF EACH CELL:

stats = regionprops3(Final_Label,'all');

CellVolumes          = stats.Volume;                   % cell volumes
CellSurfaceAreas     = stats.SurfaceArea;              % cell surface area
CellCentroids        = stats.Centroid;                 % cell centroids
CellVolSurfAreaRatio = CellVolumes./CellSurfaceAreas;       % cell volume-surface area ratio
CellEllipticity      = (stats.PrincipalAxisLength(:,1) - stats.PrincipalAxisLength(:,3))./(stats.PrincipalAxisLength(:,1)); % ellipticity
LongestAxis          = stats.PrincipalAxisLength(:,1); % length of the longest axis (diameter of longest axis)
CellElongation       = LongestAxis./((stats.PrincipalAxisLength(:,2).*stats.PrincipalAxisLength(:,3))./2); % cell elongation: longets axis divided by the average of thwo shortest
NumberOfCells        = size(stats.Volume,1);           % number of cells

%% 9.2) HISTOGRAM PLOT FOR VOLUME DISTRIBUTION
figure                                                              
hist(CellVolumes,100)
%===================================================================%
title('Nucleus volumes, number of nuclei = 526')                       
%===================================================================%
ylabel('Number of nuclei with certain volume')
xlabel('Nucleus Volume in voxels')
%% 10) Preparation for visualization of spatial parameter statistics

% 10.1) All the parameter values
stats_matrix_all      = zeros(NumberOfCells,5);                    
stats_matrix_all(:,1) = CellVolumes;                % parameter 1: volume              
stats_matrix_all(:,2) = CellVolSurfAreaRatio;       % parameter 2: vol-surfArea ratio              
stats_matrix_all(:,3) = CellEllipticity;            % parameter 3: ellipticity          
stats_matrix_all(:,4) = CellElongation;             % parameter 4: elongation
stats_matrix_all(:,5) = LongestAxis;                % parameter 5: the length of the longest axis
zscored_matrix        = zscore(stats_matrix_all);  

% 10.2) Knock out some parameters, six possibilities we may look:
stats_matrix_1thru2 = stats_matrix_all(:,1:2);        % parameters 1 and 2
stats_matrix_1thru3 = stats_matrix_all(:,1:3);        % parameters 1, 2 and 3   
stats_matrix_1thru4 = stats_matrix_all(:,1:4);        % parameters 1, 2, 3 and 4
stats_matrix_2thru5 = stats_matrix_all(:,2:5);        % parameters 2, 3, 4 and 5
stats_matrix_3thru5 = stats_matrix_all(:,3:5);        % parameters 3, 4 and 5
stats_matrix_4thru5 = stats_matrix_all(:,4:5);        % parameters 4 and 5
stats_matrix_1and3and4 = stats_matrix_all(:,[1,3,4]); % parameters 1, 3 and 4
stats_matrix_1and4 = stats_matrix_all(:,[1,4]);       % parameters 1 and 4
stats_matrix_1and3 = stats_matrix_all(:,[1,3]);       % parameters 1 and 3

% 10.3) zscore all of these partial number of spatial parameters:
zscored_1thru2    = zscore(stats_matrix_1thru2); 
zscored_1thru3    = zscore(stats_matrix_1thru3);
zscored_1thru4    = zscore(stats_matrix_1thru4);
zscored_2thru5    = zscore(stats_matrix_2thru5);
zscored_3thru5    = zscore(stats_matrix_3thru5);
zscored_4thru5    = zscore(stats_matrix_4thru5);
zscored_1and3and4 = zscore(stats_matrix_1and3and4);
zscored_1and4     = zscore(stats_matrix_1and4);
zscored_1and3     = zscore(stats_matrix_1and3);

%% 11) Save the stats-matrices into own folders for each sample: save here the label matrix and stats-matrix for the visualization purposes

save('Final_Label','Final_Label')
save('stats','stats')
save('stats_matrix_all','stats_matrix_all')

%% 12) VISUALIZATION: Clustergaram heatmaps of each single sample image: this is done once for each sample type 

% 12.1) name the parameters for heatmap in the order we have set above:
parameters           = {'Volume','Volume/Surface ratio','Ellipticity','Elongation','Longest Axis'};
parameters_1thru2    = {'Volume','Volume/Surface ratio'};
parameters_1thru3    = {'Volume','Volume/Surface ratio','Ellipticity'};
parameters_1thru4    = {'Volume','Volume/Surface ratio','Ellipticity','Elongation'};
parameters_2thru5    = {'Volume/Surface ratio','Ellipticity','Elongation','Longest Axis'};
parameters_3thru5    = {'Ellipticity','Elongation','Longest Axis'};
parameters_4thru5    = {'Elongation','Longest Axis'};
parameters_1and3and4 = {'Volume','Ellipticity','Elongation'};
parameters_1and4     = {'Volume','Elongation'};
parameters_1and3     = {'Volume','Ellipticity'};

%%
% 12.2) create heatmaps from all of these parameter-mixes
heatm_all            = clustergram(zscored_matrix','RowLabels',parameters','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_1thru2         = clustergram(zscored_1thru2','RowLabels',parameters_1thru2','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_1thru3         = clustergram(zscored_1thru3','RowLabels',parameters_1thru3','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1thru4    = clustergram(zscored_1thru4','RowLabels',parameters_1thru4','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_2thru5    = clustergram(zscored_2thru5','RowLabels',parameters_2thru5','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_3thru5    = clustergram(zscored_3thru5','RowLabels',parameters_3thru5','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_4thru5    = clustergram(zscored_4thru5','RowLabels',parameters_4thru5','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1and3and4 = clustergram(zscored_1and3and4','RowLabels',parameters_1and3and4','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1and4     = clustergram(zscored_1and4','RowLabels',parameters_1and4','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);
heatm_MEMB_1and3     = clustergram(zscored_1and3','RowLabels',parameters_1and3','ColumnPDist','cosine','RowPdist','cosine','DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);

%% 12.3) Visualization info: TELL WHAT HEATMAP DO YOU WANT TO VISUALIZE

% HERE WRITE THE CORRECT HEATMAP INFO YOU WANT TO VISUALIZE:
%===================================================================%
% 1. which heatmap? 
heatm_to_visualize = heatm_all;

% 2. what is the corresponding stat_matrix of this heatmap?
stats_matrix_to_visualize = stats_matrix_all;

% 3. which branches?
branches = [295 280]; 

% 4. what is going to be the name of the colored .tif image?
name = 'tbud-new-2branches-295-all5params.tif';
%===================================================================%
%% 12.4) Visalization itself: 
% 12.5) create an empty struct you are then filling in with parameter info needed
% to visualise the heatmap groups in original image:

xp = struct('stats_all',[],'SpatParamVals',[],'CellIdentities',[],'Centroid',[],...
    'SpatParamVals_len',[]);  

% 12.6) Fill this struct with info you told above
xp.stats_all         = load('stats');
xp.Centroid          = xp.stats_all.stats.Centroid;
xp.SpatParamVals     = stats_matrix_to_visualize;
xp.CellIdentities    = find(xp.stats_all.stats.Volume);
xp.SpatParamVals_t   = xp.SpatParamVals';
xp.SpatParamVals_len = size(xp.SpatParamVals,1);

get(heatm_to_visualize)

% 12.7) Get the heatmap you stated above and assign colormap for the
% picture: note that these are not only colours you can use.

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
cmp(10,:) = [0.9,0.4,0.5]; % baby pink       

% The counter is like an arrow pointing that from which colormap color do
% you want to start. If counter = 1, means your first colour to use is
% cmp(1,:) which is currently light blue. 
%===================================================================%
counter=1;
%===================================================================%


sub_cluster = struct('cells_of_interest',[]);
CellIdentities = cat(1,xp.CellIdentities);
range = [0;cumsum(cat(1,xp.SpatParamVals_len))];
figure
imshow(max(original_img,[],3),[]) 
hold on
for n = branches 
        group_of_interest = clusterGroup(heatm_to_visualize, n, 'col');
        Col_Labels = group_of_interest.ColumnLabels; % here we have extracted some column labels from the clustergram
        Double_Labels = cell(0);
        for i = 1 : length(Col_Labels)
            Double_Labels{i} = str2double(Col_Labels{i});
        end
        goi = cell2mat(Double_Labels);
        cells_of_Interest = CellIdentities(goi);
        sub_cluster(counter).cells_of_interest = cells_of_Interest;
    plot(xp.stats_all.stats.Centroid(sub_cluster(counter).cells_of_interest,1),...
         xp.stats_all.stats.Centroid(sub_cluster(counter).cells_of_interest,2),'*','color',cmp(counter,:),'LineWidth',3); 
    counter = counter+1;
end
hold off
sub_cluster_pruned = sub_cluster;
for i = 1 : length(sub_cluster)
    template = sub_cluster(i).cells_of_interest;
    for j = 1 : length(sub_cluster)
        pattern = sub_cluster(j).cells_of_interest;
        if isempty(setdiff(pattern,template)) &&(i~=j)
           sub_cluster_pruned(i).cells_of_interest = setdiff(template,pattern); 
        end
    end
    cells_of_Interest = sub_cluster_pruned(i).cells_of_interest; 
    in_group = intersect(find(cells_of_Interest>range(1)),find(cells_of_Interest<range(2)));      
    template = xp.CellIdentities; 
    pattern = cells_of_Interest(in_group);
    D = pdist2(template,pattern);
    in_xp_num = find(min(D,[],2)==0);
    id_group = xp.CellIdentities(in_xp_num);
    sub_cluster_pruned(i).cells_of_interest = id_group; 
end
Label_sub = 0*Final_Label;  

for c = 1 : length(sub_cluster_pruned)
    cells_of_Interest = sub_cluster_pruned(c).cells_of_interest;
    for i = 1 : length(cells_of_Interest) 
        Label_sub(Final_Label == cells_of_Interest(i)) = c;
    end
end

for z = 1 : size(Label_sub,3)
    temp  = zeros(size(Label_sub,1),size(Label_sub,2),3,'uint8');
    for c = 1 : length(sub_cluster_pruned)
        tmp = Label_sub(:,:,z);
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
