
clear all;

addpath /Users/d.deflorio/Desktop/simulation/FrequencyDiscrimination


load('spikes_young.mat','spike_y30','spk_s_y30')

n_SA1_y = 121;
n_RA_y = [122 317];
n_PC_y = [318 353];
n_aff = 353;

datafile = fopen('check_iteration', 'wt');

%% For 30 Hz
%% Compute spike train distance on individual afferent population

% - Diar - To define the temporal distribution of the spike distance metric i.e.
% tolerance for a time difference between two spikes i.e. temporal
% resolution of the afferents themselves 
cost = 500; % s^-1 = 1/cost s

ds = [spk_s_y30; spike_y30];

% compute spike distance for all afferent together and separatedly

parpool(2)
RMSyoungALL = zeros(600, 600);

dist = zeros(1,n_aff);

% Main function (compute_nomalized_dist)
parfor i = 1:size(ds,1)
    
    for j = i+1:size(ds,1)
        for aff = 1:n_aff
            
            dist(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost); % - Diar Takes for ever 
        end
        temp = dist; %temp = squeeze(temp);
        RMSyoungALL(i,j) = rssq(temp); % Diar - Outcome measure i.e. distance value for each afferent and each comparison of trials 
        fprintf(datafile,'\n%d \t', i);
        fprintf(datafile,'%d \t', j);
        
    end
end
%save('FD_RMS30Hz.mat','RMSyoungALL')
fclose (datafile);

%% for individual afferent types

RMSyoungSA1 = zeros(600, 600);

for i = 1:size(ds,1)
    for j = 1:size(ds,1)
        parfor aff = 1:n_SA1_y
            
            distSA1(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost);
            
        end
        temp = distSA1; %temp = squeeze(temp);
        RMSyoungSA1(i,j) = rms(temp);
    end
end

RMSyoungRA = zeros(600, 600);

for i = 1:size(ds,1)
    for j = 1:size(ds,1)
        parfor aff = n_RA_y(1):n_RA_y(end)
            
            distRA(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost);
            
        end
        temp =  distRA(:,n_RA_y(1):n_RA_y(end)); %temp = squeeze(temp);
        RMSyoungRA(i,j) = rms(temp);
    end
end

RMSyoungPC = zeros(600, 600);


for i = 1:size(ds,1)
    for j = 1:size(ds,1)
        parfor aff = n_PC_y(1):n_PC_y(end)
            
            distPC(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost);
            
        end
        temp = distRA(:,n_PC_y(1):n_PC_y(end)); %temp = squeeze(temp);
        RMSyoungPC(i,j) = rms(temp);
    end
end

save('FD_RMS30Hz.mat','RMSyoungALL',...
    'RMSyoungRA','RMSyoungSA1','RMSyoungPC');


%clear all;

%addpath /Users/d.deflorio/Desktop/simulation/FrequencyDiscrimination

%%
load('spikes_young.mat','spike_y300','spk_s_y300')


%% For 300 Hz
%% Compute spike train distance on individual afferent population


cost = 0; % s^-1 = 1/cost s

ds = [spk_s_y300; spike_y300];

% compute spike distance for all afferent together and separatedly

parpool(40)
RMSyoungALL = zeros(600, 600);

for i = 1:size(ds,1)
    for j = 1:size(ds,1)
        parfor aff = 1:n_aff
            
            dist(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost);
        end
        temp = dist; %temp = squeeze(temp);
        RMSyoungALL(i,j) = rms(temp);
        %fprintf(datafile,'\n%d \t', i);
        %fprintf(datafile,'%d \t', j);


    end
end
%save('FD_RMS30Hz.mat','RMSyoungALL')
%fclose (datafile);

%% for individual afferent types

RMSyoungSA1 = zeros(600, 600);

for i = 1:size(ds,1)
    for j = 1:size(ds,1)
        parfor aff = 1:n_SA1_y
            
            distSA1(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost);
            
        end
        temp = distSA1; %temp = squeeze(temp);
        RMSyoungSA1(i,j) = rms(temp);
    end
end

RMSyoungRA = zeros(600, 600);

for i = 1:size(ds,1)
    for j = 1:size(ds,1)
        parfor aff = n_RA_y(1):n_RA_y(end)
            
            distRA(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost);
            
        end
        temp =  distRA(:,n_RA_y(1):n_RA_y(end)); %temp = squeeze(temp);
        RMSyoungRA(i,j) = rms(temp);
    end
end

RMSyoungPC = zeros(600, 600);


for i = 1:size(ds,1)
    for j = 1:size(ds,1)
        parfor aff = n_PC_y(1):n_PC_y(end)
            
            distPC(aff) = compute_normalized_dist(ds{i,aff},ds{j,aff},cost);
            
        end
        temp = distRA(:,n_PC_y(1):n_PC_y(end)); %temp = squeeze(temp);
        RMSyoungPC(i,j) = rms(temp);
    end
end

save('FD_RMS300Hz.mat','RMSyoungALL',...
    'RMSyoungRA','RMSyoungSA1','RMSyoungPC');


%% The next thing is Linear Discriminant Analysis to classify these outputs 
% First split the data set before applying PCA to avoid correlations 
%  The outputs above need to be processed e.g. change in and PCA or
% multidimensional scaling the the classification can occur on the
% processed (i.e. transformed) matrix 




