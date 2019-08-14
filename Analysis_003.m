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
002_08/07/2019
 - count unstitched images
*- consider adding:
    - total image area considered
    - Image area used(contour outline of cells)
    - Mapping of lung structure
        - Bronchioles
            -peribronchiolar cells 
        - Vasculature
            -perivascular cells
        - Alveolar epithelia
003_08/12/2019
 - B: updating segmentation and adjusting thresholding
003_08/12/2019
 - B: completed segmentation and adjusting thresholding
*- B: should change fluorescence image filtering in green channel
        signal is variable within cell - perhaps aggressive gaussian
*- B: need to implement different overlap count method
        should be limited to #of cells in red channel that have red overlap
        not number of cells in green that overlap in red

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

Solution#B
    -Count the unstiched images
        -there is over-imaging: overlapping perimeter between unstitched images
%}


%% SOLUTION B
%INPUT
%{
List of images in folder
%}
%OUTPUT
%{
Data structure with data
%}

%path of images to be analyzed
%path_72hr_C = 'D:\HALI\SC8435_72h_C\Test 2_01';
file_name = '';
Im_path = '';

%path image channels 1-4
my_list = dir(Im_path);
i_increment = 1;
for i_count = 1:length(my_list)
    if (contains(my_list(i_count).name,{'CH1','CH2','CH3','CH4'}))
        Im_Data(i_increment).name = my_list(i_count).name;
        Im_Data(i_increment).folder = my_list(i_count).folder;       
        i_increment = i_increment+1;
    end
end

%select data struct
DataStruct = Im_Data;

hwait_main = waitbar(0,'Main Progress...');
hwait_main.Position = [430 460 270 56];
%process all channels individually
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
    Itemp = imadjust(Itemp);
    
    %filter
    Itemp = medfilt2(Itemp);
    Itemp_preserve = Itemp;
    
    %figure;imshow(Itemp);
    %trim outliers
    lim_high = 0.98; %eliminate 2% cells
    lim_low = 1-lim_high;
    [N_cdf_r, ~] = histcounts(Itemp,'Normalization','cdf');
    [~,lims_r(1)] = min(abs(N_cdf_r - lim_low));
    [~,lims_r(2)] = min(abs(N_cdf_r - lim_high));
    Itemp_trim = Itemp;
    Itemp_trim(Itemp_trim <= lims_r(1)) = 0;
    Itemp_trim(Itemp_trim >= lims_r(2)) = 0;
    Itemp_trim(Itemp_trim>0) = 1;
    Itemp = Itemp(logical(Itemp_trim));
    
    %set cutoff @ bottom 1/3 of hist distribution
    Icutoff = range(Itemp)/3 + min(Itemp);

    %use Icutoff level as highpass cutoff
    Itemp = Itemp_preserve;
    Itemp = Itemp(:);
    Itemp(Itemp <= Icutoff) = Icutoff;
    Itemp = reshape(Itemp,size(Itemp_preserve));
    Itemp = imadjust(Itemp);

    %binarize
    %variable fluorescence and area
    if (rem(i_count,4) == 1)
        T = adaptthresh(Itemp,0.5,'NeighborhoodSize',49,...
            'Statistic','mean');
        BW = imbinarize(Itemp,T);
    elseif rem(i_count,4) == 2
        T = adaptthresh(Itemp,0.5,'NeighborhoodSize',75,...
            'Statistic','mean');
        BW = imbinarize(Itemp,T);
    elseif rem(i_count,4) == 3
        T = adaptthresh(Itemp,0.4,'NeighborhoodSize',21,...
            'Statistic','mean');
        BW = imbinarize(Itemp,T);
    end

    %figure, imshow(BW);
    BW = bwmorph(BW,'bridge');
    BW = imfill(BW,'holes');
    BW = bwmorph(BW,'thicken',1);
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
    
    %break this iteration if no objects, must assign empty value
    if isempty(BW2_stats)
        DataStruct(i_count).Stats = BW2_stats;
    else

        %segment and stats
        BW2_stats_new = mySegmenter(BW2_stats, Itemp_preserve);
        tempCC = cell(1);
        for countObject = 1:size(BW2_stats_new)
            tempCC{1,countObject} = BW2_stats_new(countObject).PixelIdxList;
        end
        BW2_cc.NumObjects = size(BW2_stats_new,1);
        BW2_cc.PixelIdxList = tempCC;
        BW2_stats_new = regionprops(BW2_cc,{'Area','Boundingbox','Centroid','PixelIdxList'});

%         %rebuild object mask
%         BW3 = zeros(size(Itemp),'logical');
%         temp_pixelIdxList = [];
%         for k=1:length(BW2_cc.PixelIdxList)
%             temp_pixelIdxList = BW2_cc.PixelIdxList{k};
%             for j=1:length(temp_pixelIdxList)
%                 BW3(temp_pixelIdxList(j)) = 1;
%             end
%         end
%         clear tempCC
%         %make objects colorful
%         BW2_Lab = labelmatrix(BW2_cc);
%         BW2_RGB = label2rgb(BW2_Lab,'jet','k','shuffle');
%         Gridme(BW2_RGB)
% 
%         %display overlay
%         BW2_overlay = imoverlay(Itemp,BW3,mycolor);
%         Gridme(BW2_overlay)

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

% save .mat file of DataStruct
save([file_name,'_DataStruct.mat'],'DataStruct');

%count objects
num_obj = [];
for i_count = 1:length(DataStruct)   
    if rem(i_count,4) == 1
        num_obj(ceil(i_count/4),1) = length(DataStruct(i_count).Stats);        
    elseif rem(i_count,4) == 2
        num_obj(ceil(i_count/4),2) = length(DataStruct(i_count).Stats);                
    elseif rem(i_count,4) == 3
        num_obj(ceil(i_count/4),3) = length(DataStruct(i_count).Stats);
    end
end

Obj_red = sum(num_obj(:,1));
Obj_green = sum(num_obj(:,2));
Obj_blue = sum(num_obj(:,3));

%determine overlap red & green
I_Ovr = zeros(720,960,'logical');
for i_count=1:121
    r_count = (i_count*4)-3;
    g_count = (i_count*4)-2;
    Stats_tempR = DataStruct(r_count).Stats;
    Stats_tempG = DataStruct(g_count).Stats;
    %list of all objects in one image (indices)
    obj_idxR = cat(1,Stats_tempR.PixelIdxList);
    obj_idxG = cat(1,Stats_tempG.PixelIdxList);
    %overlap
    obj_idxOvr = intersect(obj_idxR,obj_idxG);
    I_Ovr(:) = 0;
    I_Ovr(obj_idxOvr) = 1;
    Ovr_CC = bwconncomp(I_Ovr);
    Ovr_obj(i_count) = Ovr_CC.NumObjects;
end
Obj_redgreen = sum(Ovr_obj);

%compute the fraction of Sftpc + Ki67 cells
% (Sftpc + Ki67 overlap) / (Sftpc)
Obj_redgreen_frac = Obj_redgreen/Obj_green;

%save cell counts (aka obj)
save([file_name,'_cellcount.mat'],...
    'Obj_red','Obj_green','Obj_blue','Obj_redgreen','Obj_redgreen_frac');
