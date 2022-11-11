clc; clear;

addpath /Users/d.deflorio/Desktop/palamedes1_10_4/Palamedes
addpath /Users/d.deflorio//Desktop/simulation/TwoPD/simulation_data
cd /Users/d.deflorio//Desktop/simulation/TwoPD/simulation_data
% % 1. Initialize parameter for classification
%

% - Diar run classification several times (rep = 50)
% rep = 50; % number of different classification
% tr_id = 200; % number of trials for each classification (out of 200)
% n_cond = 7; % number of condition
% n_tp = 15; % number of time windows

%% extract spike count
%addpath /Users/d.deflorio/Desktop/simulation/TwoPD/simulation_data
load('spikes.mat')

n_trial = 1600;
n_tw = 15; % number of windows 

n_aff = 317; 
n_aff_o = 149;
% compute distance between conditions, for each afferent
for tp = 1:n_tw
    for tr = 1:n_trial
        for aff = 1:n_aff
            spk(tp,tr,aff) =  numel(data_y{1,tp}{tr,aff});
        end
    end
end

for tp = 1:n_tw
    for tr = 1:n_trial
        for aff = 1:n_aff_o
            spk_eld(tp,tr,aff) =  numel(data_eld{1,tp}{tr,aff});
        end
    end
end

%% young
for id = 1:15
    ALL_spt{id} = squeeze(spk(id,1:200,:));
    
    ALL_2pt{id} = {squeeze(spk(id,201:400,:)); squeeze(spk(id,401:600,:)); squeeze(spk(id,601:800,:));...
        squeeze(spk(id,801:1000,:));squeeze(spk(id,1001:1200,:)); squeeze(spk(id,1201:1400,:));squeeze(spk(id,1401:1600,:))};
end

for ss = 1:7
    for tw = 1:15
        for rep = 1:50
            
            data = [ALL_spt{tw}; ALL_2pt{tw}(ss)];
            
            dt_c = [data{1}; data{2}];
            
            groups = [zeros(1,200),ones(1,200)]';
            indices = crossvalind('Kfold',groups,10);
            
            cp = classperf(groups);
            
            for i = 1:10
                test = (indices == i); train = ~test;
                
                trainingData = dt_c(train,:);
                sampleData = dt_c(test,:);
                
                Xtrain = trainingData - mean(trainingData);
                
                [coeff,score,latent,tsquared,explained,mu] = pca(trainingData);
                
                Xtest = sampleData - mu;

                tot_var = cumsum(explained);
                idx = find(tot_var > 95);
                
                if isempty(idx) == 1
                    
                    temp(i) = 0.5;
                    
                else
                    n_comp = idx(1);

                train_pca = Xtrain*coeff(:,1:n_comp);
                test_pca = Xtest*coeff(:,1:n_comp);
                
                groups = [zeros(1,180),ones(1,180)]';
                
                Mdl = fitcdiscr(train_pca, groups,'DiscrimType','pseudoLinear');
                class = Mdl.predict(test_pca);
                classperf(cp,class,test);
                
                temp(i) = cp.CorrectRate;
                end
            end
            acc_y(tw,ss,rep) = mean(temp);
        end
    end
end


%% elderly (both)
for id = 1:15
    ALL_spt{id} = squeeze(spk_eld(id,1:200,:));
    
    ALL_2pt{id} = {squeeze(spk_eld(id,201:400,:)); squeeze(spk_eld(id,401:600,:)); squeeze(spk_eld(id,601:800,:));...
        squeeze(spk_eld(id,801:1000,:));squeeze(spk_eld(id,1001:1200,:)); squeeze(spk_eld(id,1201:1400,:));squeeze(spk_eld(id,1401:1600,:))};
end

for ss = 1:7
    for tw = 1:15
        for rep = 1:50
            
            data = [ALL_spt{tw}; ALL_2pt{tw}(ss)];
            
            dt_c = [data{1}; data{2}];
            
            groups = [zeros(1,200),ones(1,200)]';
            indices = crossvalind('Kfold',groups,10);
            
            cp = classperf(groups);
            
            for i = 1:10
                test = (indices == i); train = ~test;
                
                trainingData = dt_c(train,:);
                sampleData = dt_c(test,:);
                
                Xtrain = trainingData - mean(trainingData);
                
                [coeff,score,latent,tsquared,explained,mu] = pca(trainingData);
                
                tot_var = cumsum(explained);
                idx = find(tot_var > 95);
                
                if isempty(idx) == 1
                    
                    n_comp = 1;
                    
                else
                    n_comp = idx(1);
                
                Xtest = sampleData - mu;

                train_pca = Xtrain*coeff(:,1:n_comp);
                test_pca = Xtest*coeff(:,1:n_comp);
                
                groups = [zeros(1,180),ones(1,180)]';
                
                Mdl = fitcdiscr(train_pca, groups,'DiscrimType','pseudoLinear');
                class = Mdl.predict(test_pca);
                classperf(cp,class,test);
                
                temp(i) = cp.CorrectRate;
                end
                acc_eld(tw,ss,rep) = mean(temp);

            end
        end
    end
end

%save('new_class.mat','acc_y','acc_eld')
