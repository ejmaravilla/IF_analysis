%% --------------
%UPDATES
%{
001_08/09/2019
 - validating Analysis_002 Solution B
    - effects of noise, ie artifacts
    - thresholding for masking
    - filtering of too small and too big objects
    - segmentation parameters
001_08/11/2019
 - preliminary reports
%}

%load Data
load('SC8435_72h_I_DataStruct');
DataStruct = SC8435_72h_I_DataStruct;

%find location of artifacts and validate
%manually inspect all images and select a few with artifacts
%image #31, 63, 70, 91, 92

%get list of image names
I_idx = [31, 63];
I_list = {};
count_k = 1;
for k = I_idx
    for j = 0:3
        I_ch_idx = (k*4)-3+j;
        fullfile_name = [DataStruct(I_ch_idx).folder,'\',DataStruct(I_ch_idx).name];
        fullfile_name = join(fullfile_name,'');
        I_list(count_k).ImgIdx = k;
        I_list(count_k).ImgList(j+1) = fullfile_name;
    end
    count_k = count_k+1;
end

%display raw images
for k = 1:length(I_list)
    figure('Name',['Image Index = ', num2str(I_list(k).ImgIdx)]);
    montage([I_list(k).ImgList(4),I_list(k).ImgList(1:3)], 'Size', [1 length(I_list(k).ImgList)]);
end

%display normalized images
%scale signal
I_norm = zeros(720,960,3,'uint8');
for k = 1
    for j = 1:3
        I_temp = imread(I_list(k).ImgList(j));
        %read only one channel
        if j==1
            I_temp = I_temp(:,:,1);
        elseif j==2
            I_temp = I_temp(:,:,2);
        else
            I_temp = I_temp(:,:,3);
        end
        %scale intensity range
        I_temp = imadjust(I_temp);
        %filter image
        I_temp = medfilt2(I_temp);
        %build image matrix with three channels
        if j==1
            I_norm(:,:,1) = I_temp;
        elseif j==2
            I_norm(:,:,2) = I_temp;
        else
            I_norm(:,:,3) = I_temp;
        end
    end
%    I_list(k).I_norm = I_norm;
end

%generate norm figures
temp = imread(I_list(1).ImgList(4));
montage([temp,I_norm(:,:,1),I_norm(:,:,2),I_norm(:,:,3)])

%linearize
I1temp = I_norm(:,:,1);
I2temp = I_norm(:,:,2);
I3temp = I_norm(:,:,3);
I1temp = I1temp(:);
I2temp = I2temp(:);
I3temp = I3temp(:);

%trim outliers
lim_high = 0.98;
lim_low = 1-lim_high;
[N_cdf_r, ~] = histcounts(I1temp,'Normalization','cdf');
[~,lims_r(1)] = min(abs(N_cdf_r - lim_low));
[~,lims_r(2)] = min(abs(N_cdf_r - lim_high));
[N_cdf_g, ~] = histcounts(I2temp,'Normalization','cdf');
[~,lims_g(1)] = min(abs(N_cdf_g - lim_low));
[~,lims_g(2)] = min(abs(N_cdf_g - lim_high));
[N_cdf_b, ~] = histcounts(I3temp,'Normalization','cdf');
[~,lims_b(1)] = min(abs(N_cdf_b - lim_low));
[~,lims_b(2)] = min(abs(N_cdf_b - lim_high));

I1temp_trim = I1temp;
I1temp_trim(I1temp_trim <= lims_r(1)) = 0;
I1temp_trim(I1temp_trim >= lims_r(2)) = 0;
I1temp_trim(I1temp_trim>0) = 1;
I1temp = I1temp(logical(I1temp_trim));

I2temp_trim = I2temp;
I2temp_trim(I2temp_trim <= lims_g(1)) = 0;
I2temp_trim(I2temp_trim >= lims_g(2)) = 0;
I2temp_trim(I2temp_trim>0) = 1;
I2temp = I2temp(logical(I2temp_trim));

I3temp_trim = I3temp;
I3temp_trim(I3temp_trim <= lims_b(1)) = 0;
I3temp_trim(I3temp_trim >= lims_b(2)) = 0;
I3temp_trim(I3temp_trim>0) = 1;
I3temp = I3temp(logical(I3temp_trim));

%get descriptive stats on intesity values
mean1 = mean(I1temp);
mean2 = mean(I2temp);
mean3 = mean(I3temp);
mid1 = range(I1temp)/2 + min(I1temp);
mid2 = range(I2temp)/2 + min(I2temp);
mid3 = range(I3temp)/2 + min(I3temp);
std1 = std(double(I1temp));
std2 = std(double(I2temp));
std3 = std(double(I3temp));
var1 = var(double(I1temp));
var2 = var(double(I2temp));
var3 = var(double(I3temp));

%use mid level as highpass cutoff
I1temp = I_norm(:,:,1);
I2temp = I_norm(:,:,2);
I3temp = I_norm(:,:,3);
I1temp = I1temp(:);
I2temp = I2temp(:);
I3temp = I3temp(:);

I1temp(I1temp <= floor(0.9*mid1)) = floor(0.9*mid1);
I2temp(I2temp <= mid2) = mid2;
I3temp(I3temp <= mid3) = mid3;

I1temp = reshape(I1temp,[720, 960]);
I2temp = reshape(I2temp,[720, 960]);
I3temp = reshape(I3temp,[720, 960]);
I1temp = imadjust(I1temp);
I2temp = imadjust(I2temp);
I3temp = imadjust(I3temp);

figure;imshow(I_norm(:,:,1));
figure;imshow(I1temp);
figure;imshow(I_norm(:,:,2));
figure;imshow(I2temp);
figure;imshow(I_norm(:,:,3));
figure;imshow(I3temp);

%make binary mask of objects ('cells')
BW = imbinarize(I1temp,0.15);
% figure, imshow(BW);
BW = bwmorph(BW,'bridge');
BW = imfill(BW,'holes');
BW = bwmorph(BW,'thicken',1);
BW = bwmorph(BW,'majority');
BW = bwmorph(BW,'open');
BW = bwmorph(BW,'clean');
BW1 = bwareaopen(BW,50);
figure; imshow(BW1)
% figure; imshow(I1temp)
% figure;imshowpair(I1temp,BW1,'montage')

BW = imbinarize(I2temp,0.25);
% figure, imshow(BW);
BW = bwmorph(BW,'bridge');
BW = imfill(BW,'holes');
BW = bwmorph(BW,'thicken',1);
BW = bwmorph(BW,'majority');
BW = bwmorph(BW,'open');
BW = bwmorph(BW,'clean');
BW2 = bwareaopen(BW,50);
figure; imshow(BW2)
% figure; imshow(I2temp)
% figure;imshowpair(I2temp,BW2,'montage')

BW = imbinarize(I3temp,0.35);
% figure, imshow(BW);
BW = bwmorph(BW,'bridge');
BW = imfill(BW,'holes');
BW = bwmorph(BW,'thicken',1);
BW = bwmorph(BW,'majority');
BW = bwmorph(BW,'open');
BW = bwmorph(BW,'clean');
BW3 = bwareaopen(BW,50);
figure; imshow(BW3)
% figure; imshow(I3temp)
% figure;imshowpair(I3temp,BW3,'montage')

BW1=BW3;
%connectivity map
BW1_cc = bwconncomp(BW1);
BW1_stats = regionprops(BW1_cc,{'Area','BoundingBox','Centroid','PixelIdxList'});

%segment then rebuild conncomp & stats
BW1_stats_new = mySegmenter(BW1_stats,I_norm(:,:,3));
tempCC = cell(1);
for countObject = 1:size(BW1_stats_new)
    tempCC{1,countObject} = BW1_stats_new(countObject).PixelIdxList;
end
BW1_cc.NumObjects = size(BW1_stats_new,1);
BW1_cc.PixelIdxList = tempCC;
BW1_stats_new = regionprops(BW1_cc,{'Area','Boundingbox','Centroid','PixelIdxList'});

%repeat segmentation on "large=median+3*SEM" objects
BW1_stats_new = mySegmenter_2nd(BW1_stats_new,I_norm(:,:,3));
tempCC = cell(1);
for countObject = 1:size(BW1_stats_new)
    tempCC{1,countObject} = BW1_stats_new(countObject).PixelIdxList;
end
BW1_cc.NumObjects = size(BW1_stats_new,1);
BW1_cc.PixelIdxList = tempCC;
BW1_stats_new = regionprops(BW1_cc,{'Area','Boundingbox','Centroid','PixelIdxList'});

%regenerate mask after segmentation
BW1(:) = 0;
temp_pixelIdxList = [];
for k=1:length(BW1_cc.PixelIdxList)
    temp_pixelIdxList = BW1_cc.PixelIdxList{k};
    for j=1:length(temp_pixelIdxList)
        BW1(temp_pixelIdxList(j)) = 1;
    end
end

%random color labels
BW1_Lab = labelmatrix(BW1_cc);
BW1_RGB = label2rgb(BW1_Lab,'jet','k','shuffle');
Gridme(BW1_RGB)

%single color overlay
BW1_overlay = imoverlay(imread(I_list(2).ImgList(4)),BW1,'red');
Gridme(BW1_overlay)


%generate norm figures
temp = imread(I_list(1).ImgList(4));
montage([temp,imadjust(uint8(I1_BW)),imadjust(uint8(I2_BW)),imadjust(uint8(I3_BW))])
temp = reshape([temp,temp,temp],[720 960 3]);

t = Tiff('temp.tif','w');
tagstruct.ImageLength = size(temp,1); 
tagstruct.ImageWidth = size(temp,2);
tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 3;
tagstruct.Compression = Tiff.Compression.LZW;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
tagstruct.Software = 'MATLAB';
setTag(t,tagstruct)
write(t,temp);
close(t);

t = Tiff('I1_RGB.tif','w');
tagstruct.ImageLength = size(temp,1); 
tagstruct.ImageWidth = size(temp,2);
tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 3;
tagstruct.Compression = Tiff.Compression.LZW;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
tagstruct.Software = 'MATLAB';
setTag(t,tagstruct)
write(t,I1_RGB);
close(t);

t = Tiff('I2_RGB.tif','w');
tagstruct.ImageLength = size(temp,1); 
tagstruct.ImageWidth = size(temp,2);
tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 3;
tagstruct.Compression = Tiff.Compression.LZW;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
tagstruct.Software = 'MATLAB';
setTag(t,tagstruct)
write(t,I2_RGB);
close(t);

t = Tiff('I3_RGB.tif','w');
tagstruct.ImageLength = size(temp,1); 
tagstruct.ImageWidth = size(temp,2);
tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 3;
tagstruct.Compression = Tiff.Compression.LZW;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
tagstruct.Software = 'MATLAB';
setTag(t,tagstruct)
write(t,I3_RGB);
close(t);

montage({'temp.tif','I1_RGB.tif','I2_RGB.tif','I3_RGB.tif'},'size',[1 4])



%plot Counts
figure('Name','Hist Count Red');hr = histogram(I1temp(I1temp>0)); hr.FaceColor = 'r';
hold on; xline(mean1,'--k','mean');xline(med1,'--k','median');hold off;
figure('Name','Hist Count Green');hr = histogram(I2temp); hr.FaceColor = 'g';
hold on; xline(mean2,'--k','mean');xline(med2,'--k','median');hold off;
figure('Name','Hist Count Blue');hr = histogram(I3temp); hr.FaceColor = 'b';
hold on; xline(mean3,'--k','mean');xline(med3,'--k','median');hold off;
%plot CDF
figure('Name','Hist CDF Red');hr = histogram(I1temp,'Normalization','cdf'); hr.FaceColor = 'r';
hold on; xline(mean1,'--k');hold off;
figure('Name','Hist CDF Green');hr = histogram(I2temp,'Normalization','cdf'); hr.FaceColor = 'g';
hold on; xline(mean2,'--k');hold off;
figure('Name','Hist CDF Blue');hr = histogram(I3temp,'Normalization','cdf'); hr.FaceColor = 'b';
hold on; xline(mean3,'--k');hold off;


function cdf_plot(plotdata, plotmean, plotcolor)
%make bar plot of cdf
    figure('Name',['Hist CDF ',plotcolor]);hr = histogram(plotdata,'Normalization','cdf'); hr.FaceColor = plotcolor;
    hold on; xline(plotmean,'--k');hold off;
end

function pdf_plot(plotdata, plotmean, plotcolor)
%make bar plot of pdf
    figure('Name',['Hist PDF ',plotcolor]);hr = histogram(plotdata,'Normalization','pdf'); hr.FaceColor = plotcolor;
    hold on; xline(plotmean,'--k');hold off;
end

function BW = binarize_me(Im,level)
%generate binary mask
    BW = imbinarize(Im,level);
    % figure, imshow(BW);
    BW = bwmorph(BW,'bridge');
    BW = imfill(BW,'holes');
    BW = bwmorph(BW,'thicken',1);
    BW = bwmorph(BW,'majority');
    BW = bwmorph(BW,'open');
    BW = bwmorph(BW,'clean');
    BW = bwareaopen(BW,50);
end

