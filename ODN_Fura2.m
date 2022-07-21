%% Clean up
clear;
close all
clc
% lp = 0;
% imgshow = 1;
Save = 0;
%gcp;                                            % start a parallel engine

%% Load files
[filename,foldername] = uigetfile('*.zip','Select a ROI.zip file'); % Load the txt file
strROIArchiveFilename = strcat(foldername,filename);
[cvsROIs] = ReadImageJROI(strROIArchiveFilename);

e = strfind(filename,' e');
embryo_no = strcat('0',filename(e+2));
h = strfind(filename,'_h');
hemi_no = strcat('0',filename(h+2));
exp = strfind(filename,'_0');
exp_no = strcat('0',(filename(exp+2)));

cd(strcat(foldername,'\',embryo_no));
tif_name = strcat(hemi_no,'_',exp_no,'.tif');
data = read_file(tif_name);
FOV = [size(data,1) size(data,2)];
NOF = size(data,3);

%% Define the parameters 
Bin_s = 1;   % Spatial binning factor
Bin_t = 2;   % Temporal binning factor

FOV_s = FOV/Bin_s;
numplanes = NOF/Bin_t;

N = 2;                               %   First N-point to be removed due to acquisition artifacts
P = 4;                               %   P-point after puffing to be removed due to artifacts
puff_onset = [100/Bin_t 200/Bin_t 300/Bin_t];    %   # of frame

SR = 10;                             %   Hz
SR_t = 10/Bin_t;                     %   Sampling rate after binning
dt = 1/SR_t;                         %   second

time = dt:dt:numplanes/SR_t;
time = time-puff_onset(1)/SR_t;      %  time relative to puff
time_base = dt*(N/Bin_t+1):dt:puff_onset(1)/SR_t;
time_base = time_base-puff_onset(1)/SR_t; 
time_postpuff = dt*(puff_onset(1)+P/Bin_t+1):dt:puff_onset(2)/SR_t;
time_postpuff = time_postpuff-puff_onset(1)/SR_t; 
time_c = horzcat(time_base,time_postpuff);

% fs = SR_t;
% fpass = 2;                           % lowpass filter cutoff frequency (ex: 2 Hz)
%% Spatial and/or temporal averaging
if Bin_s > 1 || Bin_t > 1
    data_avg = imresize3(data,[FOV_s numplanes]);
else
    data_avg = data;
end

%% Extract raw signals from ROIs
Data_avg_stretch = reshape(data_avg, FOV_s(1)*FOV_s(2), []);
Data_avg_stretch_bg = reshape(data_avg(:,1:0.5*FOV_s(2),:), FOV_s(1)*FOV_s(2)*0.5, []);
Data_avg_stretch_bg_sort = sort(Data_avg_stretch_bg,1);
bp = fix(FOV_s(1)*FOV_s(2)*0.5/100); % number of background pixels
ROIs.background = mean(Data_avg_stretch_bg_sort(1:bp,:));
data_temp = data_avg(:,:,1); 

figure ('Name', 'ROIs');
imagesc(data_temp)
colormap (gray)
hold on   
for n = 1:size(cvsROIs,2)   
    ROIs.bounds{n,1} = cvsROIs{1,n}.vnRectBounds/Bin_s;
    ROIs.center{n,1} = [ROIs.bounds{n,1}(2)+(ROIs.bounds{n,1}(4)-ROIs.bounds{n,1}(2))/2 ROIs.bounds{n,1}(1)+(ROIs.bounds{n,1}(3)-ROIs.bounds{n,1}(1))/2];  
    ROIs.ellipse{n,1} = images.roi.Ellipse(gca,'Center',ROIs.center{n,1},'Semiaxes',[ROIs.bounds{n,1}(4)-ROIs.center{n,1}(1) ROIs.bounds{n,1}(3)-ROIs.center{n,1}(2)]);                
    ROImask = createMask(ROIs.ellipse{n,1});
    mfROIPixels = Data_avg_stretch(reshape(ROImask, FOV_s(1)*FOV_s(2), []), :);
    ROIs.data{n,1} = mean(mfROIPixels);    
    if n == size(cvsROIs,2)
       ROIs.data{n+1,1} = ROIs.background;        
%         ROIs.background{1} = images.roi.Rectangle(gca,'Position',[0,FOV_s(1)*0.2,FOV_s(2)*0.2,FOV_s(1)*0.6],'Color','r');   % [xmin, ymin, width, height] xmin ymin: upper left corner 
%         ROImask = createMask(ROIs.background{1});
%         mfROIPixels = Data_avg_stretch(reshape(ROImask, FOV_s(1)*FOV_s(2), []), :);
%         ROIs.data{n+1,1} = mean(mfROIPixels);       
    end  
end

figure ('Name', 'Raw traces');
xlabel('Time relative to puff(s)');
hold on
for i = 1:size(cvsROIs,2)
    plot(time,ROIs.data{i})       
    if i == size(cvsROIs,2)
        plot(time,ROIs.data{i+1},'k','LineWidth',2)
    end
end

%% Signal processing
% if lp                           % low-pass filtering
%     figure ('Name', 'Low-pass filtered');
%     hold on
%     for j = 1:size(cvsROIs,2)
%         ROIs.data_lp{j,1} = lowpass(ROIs.data{j},fpass,fs);
%         plot(ROIs.data_lp{j,1})
%         if j == size(cvsROIs,2)
%             ROIs.data_lp{j+1,1} = lowpass(ROIs.data{j+1,1},fpass,fs);
%             plot(ROIs.data_lp{j+1,1},'k','LineWidth',2)
%         end
%     end
% end
                            
for i = 1:size(ROIs.data,1)
    ROIs.data_base{i,1} = ROIs.data{i}(1+N/Bin_t:puff_onset(1));      % removal of first N-point (baseline)
    ROIs.data_postpuff{i,1} = ROIs.data{i}(puff_onset(1)+P/Bin_t+1:puff_onset(2));    % removal of first P-point (postpuff)
end
                     
repeats = 20;                                         %  moving average repeats
movavg = 5;                                           %  local k-point mean
figure ('Name', 'Moving averaged');
xlabel('Time relative to puff(s)');
hold on
for k = 1:size(cvsROIs,2)
    i = 1;
    while i < repeats
        ROIs.data_base_ma{k,1} = movmean(ROIs.data_base{k},movavg);
        ROIs.data_postpuff_ma{k,1} = movmean(ROIs.data_postpuff{k},movavg);
        i = i+1;
    end
    ROIs.data_update{k,1} = horzcat(NaN(1,N/Bin_t),ROIs.data_base_ma{k,1},NaN(1,P/Bin_t),ROIs.data_postpuff_ma{k,1});
    plot(time(1:puff_onset(2)),ROIs.data_update{k,1})
    if k == size(cvsROIs,2)
        i = 1;
        while i < repeats
            ROIs.data_base_ma{k+1,1} = movmean(ROIs.data_base{k+1},movavg);
            ROIs.data_postpuff_ma{k+1,1} = movmean(ROIs.data_postpuff{k+1},movavg);
            i = i+1;
        end
        ROIs.data_update{k+1,1} = horzcat(NaN(1,N/Bin_t),ROIs.data_base_ma{k+1,1},NaN(1,P/Bin_t),ROIs.data_postpuff_ma{k+1,1});
        plot(time(1:puff_onset(2)),ROIs.data_update{k+1,1},'k','LineWidth',2)
    end
end
 
figure ('Name', 'BG_subtracted');
xlabel('Time relative to puff(s)');
hold on
for k = 1:size(cvsROIs,2)
    ROIs.data_p{k,1} = ROIs.data_update{k}-ROIs.data_update{end};      % background subtraction
    plot(time(1:puff_onset(2)),ROIs.data_p{k,1})
end

%% dF/F0 calculation
for k = 1:size(cvsROIs,2)
    ROIs.data_base_ma_mean{k,1} = mean(ROIs.data_p{k,1}(1+N/Bin_t:puff_onset(1)));          % average the baseline 
end
figure ('Name', 'dF/F0');
xlabel('Time relative to puff(s)');
hold on
for k = 1:size(cvsROIs,2)
    ROIs.data_p{k,2} = ROIs.data_p{k,1}/ROIs.data_base_ma_mean{k,1};      % normalize the baseline
    plot(time(1:puff_onset(2)),ROIs.data_p{k,2})
end

%% Linear fitting and subtraction of baseline trend
f1 = figure ('Name', 'Baseline_subtracted');
x0 = 100;
y0 = 250;
width = 800;
height = 800;
f1.Position = [x0,y0,width,height]; 

for k = 1:size(cvsROIs,2)    
    [fit.pf{k,1},fit.S{k,1}] = polyfit(time_base,ROIs.data_p{k,2}(1+N/Bin_t:puff_onset(1)),1);   % 1+N/Bin_t
    [fit.y_fit{k,1},fit.delta{k,1}] = polyval(fit.pf{k},time(1:puff_onset(2)),fit.S{k});
    ROIs.data_p{k,3} = ROIs.data_p{k,2}-fit.y_fit{k,1};
    subplot(2,1,1)
    plot(time(1:puff_onset(2)),ROIs.data_p{k,2})
    plot(time(1:puff_onset(2)),fit.y_fit{k,1},'r')
    title('Baseline fitting')
    xticklabels({});
    hold on
    
    subplot(2,1,2)
    plot(time(1:puff_onset(2)),ROIs.data_p{k,3})
    title('Baseline subtracted')
    xlabel('Time relative to puff(s)');
    hold on
%     plot(time_base,fit.y_fit{k,1}+2*fit.delta{k,1},'m--',time_base,fit.y_fit{k,1}-2*fit.delta{k,1},'m--')    % 95% Prediction Interval     
end

%% Quality controls, ROI exclusion
SD_factor = 3;
Onset_thre = 10;
Onset_avg = 5;
Onset_avg = Onset_avg-1;
Onset_ind = puff_onset(1)+P/Bin_t+1;
f2 = figure ('Name', 'ROI exclusion');
f2.Position = [x0,y0,width,height]; 

for k = 1:size(cvsROIs,2)
    ROIs.base_SD{k,1} = std(ROIs.data_p{k,3}(1+N/Bin_t:puff_onset(1)));
    ROIs.postpuff_median{k,1} = median(ROIs.data_p{k,3}(Onset_ind:puff_onset(2)));
    
    if abs(mean(ROIs.data_p{k,3}(Onset_ind:Onset_ind+Onset_avg))) > Onset_thre*ROIs.base_SD{k,1} && ROIs.postpuff_median{k,1} < mean(ROIs.data_p{k,3}(Onset_ind:Onset_ind+Onset_avg))  ...
         || mean(ROIs.data_p{k,3}(Onset_ind:Onset_ind+Onset_avg)) < -Onset_thre*ROIs.base_SD{k,1}
        subplot(3,1,3)
        plot(time(1:puff_onset(2)),ROIs.data_p{k,3})
        title('Excluded ROIs')
        xlabel('Time relative to puff(s)');
        hold on
    else
        if ROIs.postpuff_median{k,1} < SD_factor*ROIs.base_SD{k,1}
            ROIs.ROIs_NR{k,1} = ROIs.data_p{k,3};
            subplot(3,1,2)
            plot(time(1:puff_onset(2)),ROIs.data_p{k,3})
            title('Non-responding ROIs')
            xlabel('Time relative to puff(s)');
            hold on
        else 
            ROIs.ROIs_R{k,1} = ROIs.data_p{k,3};
            subplot(3,1,1)
            plot(time(1:puff_onset(2)),ROIs.data_p{k,3})
            title('Responding ROIs')
            xlabel('Time relative to puff(s)');
            hold on
        end
    end
end
%% data reformatting and mean calculation
line_wd = 1.5;
if isfield(ROIs, 'ROIs_NR')   
    ROIs.ROIs_NR_merge = vertcat(ROIs.ROIs_NR{:});
    ROIs.ROIs_NR_merge_mean = mean(ROIs.ROIs_NR_merge);
    ROIs.numROI_NR = size(ROIs.ROIs_NR_merge,1);
    subplot(3,1,2)
    plot(time(1:puff_onset(2)),ROIs.ROIs_NR_merge_mean,'k','LineWidth',line_wd)
end
if isfield(ROIs, 'ROIs_R')    
    ROIs.ROIs_R_merge = vertcat(ROIs.ROIs_R{:});
    ROIs.ROIs_R_merge_mean = mean(ROIs.ROIs_R_merge);
    ROIs.numROI_R = size(ROIs.ROIs_R_merge,1);
    subplot(3,1,1)
    plot(time(1:puff_onset(2)),ROIs.ROIs_R_merge_mean,'k','LineWidth',line_wd)
end

%% Area under curve calculation
if isfield(ROIs, 'ROIs_R')  
    ROIs.ROIs_R_merge_filtered = ROIs.ROIs_R_merge;
    ROIs.base_SD_ind = find(~cellfun(@isempty,ROIs.ROIs_R));
    for k = 1:ROIs.numROI_R
        AUC_ind = ROIs.ROIs_R_merge(k,Onset_ind:puff_onset(2)) < SD_factor*ROIs.base_SD{ROIs.base_SD_ind(k),1};
        AUC_ind2(Onset_ind:puff_onset(2)) = AUC_ind;
        ROIs.ROIs_R_merge_filtered(k,AUC_ind2) = 0;
        ROIs.ROIs_R_AUC(k,1) = trapz(ROIs.ROIs_R_merge_filtered(k,Onset_ind:puff_onset(2)));
    end
    if ~isnan(ROIs.ROIs_R_merge_mean) | ~isempty(ROIs.ROIs_R_merge_mean)
        ROIs.base_mean_SD = std(ROIs.ROIs_R_merge_mean(1,1+N/Bin_t:puff_onset(1)));
        ROIs.ROIs_R_merge_mean(2,:) = ROIs.ROIs_R_merge_mean(1,:);
        AUC_ind = ROIs.ROIs_R_merge_mean(1,Onset_ind:puff_onset(2)) < SD_factor*ROIs.base_mean_SD;
        AUC_ind2(Onset_ind:puff_onset(2)) = AUC_ind;
        ROIs.ROIs_R_merge_mean(2,AUC_ind2) = 0;
        ROIs.ROIs_R_mean_AUC = trapz(ROIs.ROIs_R_merge_mean(2,Onset_ind:puff_onset(2)));
    end
end

%% data exporting
if Save
    jpg_name = strcat(hemi_no,'_',exp_no,'.jpeg');
    saveas(f2,jpg_name)

    if isfield(ROIs, 'ROIs_R')  
        txt_name = strcat(hemi_no,'_',exp_no,'_ROIs_R.txt');
        writematrix(ROIs.ROIs_R_merge', txt_name)
    end

    if isfield(ROIs, 'ROIs_NR')  
        txt_name = strcat(hemi_no,'_',exp_no,'_ROIs_NR.txt');
        writematrix(ROIs.ROIs_NR_merge', txt_name)
    end
end

if isfield(ROIs, 'ROIs_R')  
    n = ROIs.numROI_R;
else
    n = 1;
end
sz = [n 5];
varTypes = {'doublenan','doublenan','doublenan','doublenan','doublenan'};
varNames = {'#ROI_NR' '#ROI_R' '%ROI_R' 'ROI_R_AUC (%*s)' 'ROI_R_mean_AUC (%*s)'};
T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

if isfield(ROIs, 'ROIs_NR')  
    T{1,'#ROI_NR'} = ROIs.numROI_NR;
    if isfield(ROIs, 'ROIs_R')  
        T{1,'%ROI_R'} = 100*ROIs.numROI_R/(ROIs.numROI_NR+ROIs.numROI_R);
    else
        T{1,'%ROI_R'} = 0; 
    end
else
    if isfield(ROIs, 'ROIs_R')  
        T{1,'%ROI_R'} = 100;
    else
        T{1,'%ROI_R'} = NaN;
    end
end

if isfield(ROIs, 'ROIs_R')  
    T{1,'#ROI_R'} = ROIs.numROI_R;
    T{:,'ROI_R_AUC (%*s)'} = ROIs.ROIs_R_AUC;
    if ~isnan(ROIs.ROIs_R_merge_mean) | ~isempty(ROIs.ROIs_R_merge_mean)
        T{1,'ROI_R_mean_AUC (%*s)'} = ROIs.ROIs_R_mean_AUC;
    end
end
if Save
    xlsx_name = strcat(hemi_no,'_',exp_no,'.xlsx');
    writetable(T,xlsx_name);
end


