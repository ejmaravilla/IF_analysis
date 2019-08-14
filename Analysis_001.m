%% --------------
%UPDATES
%{
001_08/06/2019
 - setup and test analysis
*- need to implement photoshop preprocessing in matlab
001_08/06/2019
*- need to work on analyzing individual images and computing #overlap
*- need to segment all cells in view using dapi and brightfield
*- need to workout stitching individual channels
%}

%% --------------
%Objective
%{
Analyze IHC fluorescence images
Images captures on Keyence Fluorescent Microscope
    40x
    Files:
        - stitched RGB image of field
        - unstitched unordered images with three channels
            -R = Ki67
            -G = Sftpc
            -B = dapi
Treatment Conditions
    -0hr hyperoxia
    -72hr hyperoxia

Solution#A
    -Use stitched RGB image of entire field to count:
        %(Ki67+Sftpc) = (#Ki67 + #Sftpc)/#Ki67

Solution#B
    -Count the unstiched images
        -there is over-imaging: overlapping perimeter between unstitched images

Solution#C
    -Determine sorting pattern and overlap correction to 
        stitch individual channels together & count cells
        in each respective channel
%}

%% SOLUTION A
%INPUT
%{
     *RGB image
    Greyscale RGB image with cells of interest have very dim pixel
    intensities
%}
%OUTPUT
%{
    csv table
%}

%add path
    %done manually
%load images
file_a1 = 'Ki67_0hr_152333.tif';
file_a2 = 'Sftpc_0hr_152333.tif';
file_a3 = 'combo_Ki67_sftpc_0hr_152333.tif';
file_e1 = 'Ki67_72hr_143359.tif';
file_e2 = 'Sftpc_72hr_143359.tif';
file_e3 = 'combo_Ki67_Sftpc_72hr_143359.tif';
%imshow(control_ki67_0hr)

%file list
file_list = {'Ki67_0hr_152333.tif';
'Sftpc_0hr_152333.tif';
'combo_Ki67_sftpc_0hr_152333.tif';
'Ki67_72hr_143359.tif';
'Sftpc_72hr_143359.tif';
'combo_Ki67_Sftpc_72hr_143359.tif'};

%loop through each file
for i = 1:length(file_list)
    ImageData(i).Name = file_list{i};
    I = imread(ImageData(i).Name);
    %convert to greyscale single channel
    I_temp = zeros(size(I(:,:,1)),'uint8');
    for j=1:3
        I_temp = I_temp + I(:,:,j);
    end

    I_temp2 = rgb2gray(I);

    %scale image
    I_temp = double(I_temp);
    I_temp = I_temp - min(I_temp(:));
    I_temp = uint8(I_temp*256/(max(I_temp(:))));

    I_temp2 = double(I_temp2);
    I_temp2 = I_temp2 - min(I_temp2(:));
    I_temp2 = uint8(I_temp2*256/(max(I_temp2(:))));

    %trim border values by 5 pixels
    I_temp(1:10,:) = 256;
    I_temp(:,1:10) = 256;
    I_temp(5665:5674,:) = 256;
    I_temp(:,7565:7571) = 256;

    I_temp2(1:10,:) = 256;
    I_temp2(:,1:10) = 256;
    I_temp2(5665:5674,:) = 256;
    I_temp2(:,7565:7571) = 256;

    %binary
    BW = ~im2bw(I_temp,0.55);
    %figure, imshow(BW);
    BW = bwmorph(BW,'bridge');
    BW = imfill(BW,'holes');
    BW = bwmorph(BW,'clean');
    BW = bwareaopen(BW,16);
    BW2 = BW;
    %figure; imshow(BW2)

    %morph
    BW2 = bwmorph(BW2,'close');
    BW2 = bwareaopen(BW2,25);
    %figure;imshowpair(I_temp,BW2,'montage')

    %connectivity map
    BW2_cc = bwconncomp(BW2);
    BW2_stats = regionprops(BW2_cc,{'Area','BoundingBox','Centroid','PixelIdxList'});

    I_temp_inverse = imadjust(I_temp2,[0;1],[1;0]); %invert image
    BW2_stats_new = mySegmenter(BW2_stats,I_temp_inverse);

    %dislply label matrix after watershed segmentation
    tempCC = cell(1);
    for countObject = 1:size(BW2_stats_new)
        tempCC{1,countObject} = BW2_stats_new(countObject).PixelIdxList;
    end
    BW2_cc.NumObjects = size(BW2_stats_new,1);
    BW2_cc.PixelIdxList = tempCC;
    BW2_stats_new = regionprops(BW2_cc,{'Area','Boundingbox','Centroid','PixelIdxList'});

    BW3 = zeros(size(I_temp),'logical');
    temp_pixelIdxList = [];
    for k=1:length(tempCC)
        temp_pixelIdxList = tempCC{k};
        for j=1:length(temp_pixelIdxList)
            BW3(temp_pixelIdxList(j)) = 1;
        end
    end
    clear tempCC

    BW2_Lab = labelmatrix(BW2_cc);
    BW2_RGB = label2rgb(BW2_Lab,'jet','k','shuffle');
%     Gridme(BW2_RGB)

    %display overlay after watershed
    BW2_overlay = imoverlay(I_temp2,BW3,'red');
%     Gridme(BW2_overlay)
    % another way to get an overlay image
    % I5adjust_red = I5adjust(:,:,1);
    % for k = 1:length(I5_ki67_stats)
    %         I5adjust_red(I5_ki67_stats(k).PixelIdxList) = 255;
    % end
    % I5adjust_ki67_overlay = I5adjust;
    % I5adjust_ki67_overlay(:,:,1) = I5adjust_red.*uint8(MaskCombo5);
    % Gridme(I5adjust_ki67_overlay)
    
    ImageData(i).Stats = BW2_stats_new;
    ImageData(i).I_inverse = I_temp_inverse;
    ImageData(i).I_RGB = BW2_RGB;
    ImageData(i).I_overlay = BW2_overlay;
    
end

%generate segmentation on all cells
%map fluorescence signals to cells
%determine cell overlap by ki67 and sftpc signals

%calculations
perc_Ki67_0hr = (5900)/862 * 100;
perc_Ki67_72hr = (11434)/5858 * 100;

%validation
combo_72hr = ImageData(6).I_RGB;
combo_0hr = ImageData(3).I_RGB;

Gridme(combo_0hr)
Gridme(combo_72hr)


%% SOLUTION B
%INPUT
%{
List of images in folder
%}
%OUTPUT
%{
Data structure with data
%}

%path of 72hr treatment
path_72hr = 'D:\HALI\SC8435_72h_I\Test1';
path_0hr = 'D:\HALI\SC8507_0h_C\Test1_01';

%path_72hr
my_list = dir(path_72hr);
i_increment = 1;
for i_count = 1:length(my_list)
    if (contains(my_list(i_count).name,{'CH1','CH2','CH3','CH4'}))
        SC8435_72h_I_Data(i_increment).name = my_list(i_count).name;
        SC8435_72h_I_Data(i_increment).folder = my_list(i_count).folder;       
        i_increment = i_increment+1;
    end
end

%path_0hr
my_list = dir(path_0hr);
i_increment = 1;
for i_count = 1:length(my_list)
    if (contains(my_list(i_count).name,{'CH1','CH2','CH3','CH4'}))
        SC8507_0h_C_Data(i_increment).name = my_list(i_count).name;
        SC8507_0h_C_Data(i_increment).folder = my_list(i_count).folder;
        i_increment = i_increment+1;
    end
end

%select data struct
DataStruct = SC8435_72h_I_Data;
DataStruct = SC8507_0h_C_Data;

hwait_main = waitbar(0,'Main Progress...');
hwait_main.Position = [430 460 270 56];
%process all channels as binary images
for i_count = 1:length(DataStruct)
    fileName = [DataStruct(i_count).folder,'\',DataStruct(i_count).name];
    Itemp = imread(fileName);
    %keep one channel for fluorescence
    if (rem(i_count,4) == 1)
        Itemp = Itemp(:,:,1); %red ki67
        mycolor = 'red';
    elseif rem(i_count,4) == 2
        Itemp = Itemp(:,:,2); %green Sftpc
        mycolor = 'green';
    elseif rem(i_count,4) == 3
        Itemp = Itemp(:,:,3); %blue dapi
        mycolor = 'blue';
    end
    
    %scale signal
    Itemp = double(Itemp);
    Itemp = Itemp - min(Itemp(:));
    Itemp = uint8(Itemp*256/(max(Itemp(:))));
    
    %filter
    Itemp = medfilt2(Itemp);
%     figure;imshow(Itemp);
    
    %binarize
    BW = imbinarize(Itemp,0.2);
    %figure, imshow(BW);
    BW = bwmorph(BW,'bridge');
    BW = imfill(BW,'holes');
    BW = bwmorph(BW,'thicken',3);
    BW = bwmorph(BW,'majority');
    BW = bwmorph(BW,'open');
    BW = bwmorph(BW,'clean');
    BW = bwareaopen(BW,50);
    BW2 = BW;
%     figure; imshow(BW2)
%     figure; imshow(Itemp)
%     figure;imshowpair(Itemp,BW2,'montage')

    %connectivity map
    BW2_cc = bwconncomp(BW2);
    BW2_stats = regionprops(BW2_cc,{'Area','BoundingBox','Centroid','PixelIdxList'});
    
    %break this iteration if no objects
    if isempty(BW2_stats)
        DataStruct(i_count).Stats = BW2_stats;
    else

        %get the brightfield image for segmentation
        bf_idx = ceil(i_count/4)*4; %every 4th image in list
        I_bf = imread(DataStruct(bf_idx).name);

        %segment and stats
        I_bf_inverse = imadjust(I_bf,[0;1],[1;0]); %invert image
        BW2_stats_new = mySegmenter(BW2_stats,I_bf_inverse);

        %generate binary image
        %dislply label matrix after watershed segmentation
        tempCC = cell(1);
        for countObject = 1:size(BW2_stats_new)
            tempCC{1,countObject} = BW2_stats_new(countObject).PixelIdxList;
        end
        BW2_cc.NumObjects = size(BW2_stats_new,1);
        BW2_cc.PixelIdxList = tempCC;
        BW2_stats_new = regionprops(BW2_cc,{'Area','Boundingbox','Centroid','PixelIdxList'});

    %     BW3 = zeros(size(Itemp),'logical');
    %     temp_pixelIdxList = [];
    %     for k=1:length(tempCC)
    %         temp_pixelIdxList = tempCC{k};
    %         for j=1:length(temp_pixelIdxList)
    %             BW3(temp_pixelIdxList(j)) = 1;
    %         end
    %     end
    %     clear tempCC

    %     BW2_Lab = labelmatrix(BW2_cc);
    %     BW2_RGB = label2rgb(BW2_Lab,'jet','k','shuffle');
    %     Gridme(BW2_RGB)

        %display overlay after watershed
    %     BW2_overlay = imoverlay(I_bf,BW3,mycolor);
    %     Gridme(BW2_overlay)

        %assign to data structure
        DataStruct(i_count).Stats = BW2_stats_new;
    %     DataStruct(i_count).Im_BW = BW3;
    %     DataStruct(i_count).Im_RGB = BW2_RGB;
    %     DataStruct(i_count).Im_overlay = BW2_overlay;
    end
    %progress bar
    hwait_main.Color = [mod(i_count,2) 1 mod(i_count,2)];
    waitbar(i_count/length(DataStruct),hwait_main)
    
end
delete(hwait_main)

%redesignate names
% SC8507_0h_C_DataStruct = DataStruct;
% SC8435_72h_I_DataStruct = DataStruct;

num_obj = [];
for i_count = 1:length(SC8507_0h_C_DataStruct)   
    if rem(i_count,4) == 1
        num_obj(ceil(i_count/4),1) = length(SC8507_0h_C_DataStruct(i_count).Stats);        
    elseif rem(i_count,4) == 2
        num_obj(ceil(i_count/4),2) = length(SC8507_0h_C_DataStruct(i_count).Stats);                
    elseif rem(i_count,4) == 3
        num_obj(ceil(i_count/4),3) = length(SC8507_0h_C_DataStruct(i_count).Stats);
    end
end

SC8507_0h_C_red = sum(num_obj(:,1));
SC8507_0h_C_green = sum(num_obj(:,2));
SC8507_0h_C_blue = sum(num_obj(:,3));

num_obj = [];
for i_count = 1:length(SC8435_72h_I_DataStruct)   
    if rem(i_count,4) == 1
        num_obj(ceil(i_count/4),1) = length(SC8435_72h_I_DataStruct(i_count).Stats);        
    elseif rem(i_count,4) == 2
        num_obj(ceil(i_count/4),2) = length(SC8435_72h_I_DataStruct(i_count).Stats);                
    elseif rem(i_count,4) == 3
        num_obj(ceil(i_count/4),3) = length(SC8435_72h_I_DataStruct(i_count).Stats);
    end
end

SC8435_72h_I_red = sum(num_obj(:,1));
SC8435_72h_I_green = sum(num_obj(:,2));
SC8435_72h_I_blue = sum(num_obj(:,3));

%determine overlap red & green
DataStruct=SC8507_0h_C_DataStruct;
I_Ovr = zeros(720,960,'logical');
for i_count=1:121
    r_count = (i_count*4)-3;
    g_count = (i_count*4)-2;
    Stats_tempR = DataStruct(r_count).Stats;
    Stats_tempG = DataStruct(g_count).Stats;
    %list of indices
    obj_idxR = cat(1,Stats_tempR.PixelIdxList);
    obj_idxG = cat(1,Stats_tempG.PixelIdxList);
    %overlap
    obj_idxOvr = intersect(obj_idxR,obj_idxG);
    I_Ovr(:) = 0;
    I_Ovr(obj_idxOvr) = 1;
    Ovr_CC = bwconncomp(I_Ovr);
    Ovr_obj(i_count) = Ovr_CC.NumObjects;
end
SC8507_0h_C_redgreen = sum(Ovr_obj);

DataStruct=SC8435_72h_I_DataStruct;
I_Ovr = zeros(720,960,'logical');
for i_count=1:121
    r_count = (i_count*4)-3;
    g_count = (i_count*4)-2;
    Stats_tempR = DataStruct(r_count).Stats;
    Stats_tempG = DataStruct(g_count).Stats;
    %list of indices
    obj_idxR = cat(1,Stats_tempR.PixelIdxList);
    obj_idxG = cat(1,Stats_tempG.PixelIdxList);
    %overlap
    obj_idxOvr = intersect(obj_idxR,obj_idxG);
    I_Ovr(:) = 0;
    I_Ovr(obj_idxOvr) = 1;
    Ovr_CC = bwconncomp(I_Ovr);
    Ovr_obj(i_count) = Ovr_CC.NumObjects;
end
SC8435_72h_I_redgreen = sum(Ovr_obj);


%% SOLUTION C
%INPUT
%{

%}
%OUTPUT
%{

%}


%% MISC

%texture segmenting
I_texture_seg = entropyfilt(I_temp2);
Eim = mat2gray(I_texture_seg);
imshow(Eim);
BW_texture_seg = imbinarize(Eim, .8);
imshow(BW_texture_seg);
