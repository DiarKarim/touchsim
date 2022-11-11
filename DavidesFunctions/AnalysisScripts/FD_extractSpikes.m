%% extract spike for young
%% 300 Hz

cd /Users/d.deflorio/Desktop/simulation/FrequencyDiscrimination

load('afferent')

n_SA1_y = sum(aff_y.iSA1);
n_RA_y = sum(aff_y.iRA);
n_PC_y = sum(aff_y.iPC);
n_aff_y = n_SA1_y + n_RA_y + n_PC_y;

n_SA1_old = sum(aff_old.iSA1);
n_RA_old = sum(aff_old.iRA);
n_PC_old = sum(aff_old.iPC);
n_aff_old = n_SA1_old + n_RA_old + n_PC_old;

load('FD_300Hz_young.mat');

% Extract spike trains
n_sim = 150;
cond = 5;
n_aff = n_aff_y;

r_standard_y300 = reshape(r_standard_y300, 1, n_sim);
r_comp_y300 = reshape(permute(r_comp_y300, [1 3 2]), cond, n_sim);

% plot(r_comp_y3002{1,101}.stimulus.trace); hold on;
% plot(r_comp_y3002{5,101}.stimulus.trace,'r');
% plot(r_comp_y3002{1,1}.stimulus.trace,'g');

% young
for tt = 1:cond
    for ss = 1:n_sim
        for aff = 1:n_aff
            
            spk_y300{tt,ss,aff} = r_comp_y300{tt,ss}.responses(1,aff).spikes;
            spk_s_y300{ss,aff} =  r_standard_y300{ss}.responses(1,aff).spikes;
            
        end
    end
end
% plot(r_standard_y300{1,1}.stimulus); 
% figure(2)
% plot(r_comp_y3002{3,1}.stimulus);

% reshape spike matrices
spike_y300 = reshape(permute(spk_y300,[2 1 3]), cond*n_sim, n_aff);

clear r_standard_y300 r_comp_y300

% 30 Hz
load('FD_30Hz_young.mat');

% Extract spike trains
n_sim = 150;
cond = 5;

r_standard_y30 = reshape(r_standard_y30, 1, n_sim);
r_comp_y30 = reshape(permute(r_comp_y30, [1 3 2]), cond, n_sim);


% young
for tt = 1:cond
    for ss = 1:n_sim
        for aff = 1:n_aff
            
            spk_y30{tt,ss,aff} = r_comp_y30{tt,ss}.responses(1,aff).spikes;
            spk_s_y30{ss,aff} =  r_standard_y30{ss}.responses(1,aff).spikes;
            
        end
    end
end


% reshape spike matrices
spike_y30 = reshape(permute(spk_y300,[2 1 3]), cond*n_sim, n_aff);


save('spikes_young.mat','spk_s_y300','spike_y300','spk_s_y30','spike_y30')

clear r_standard_y30 r_comp_y30

%% Elderly 300 Hz
load('FD_300Hz_old.mat');

% Extract spike trains
n_sim = 150;
cond = 5;

r_standard_old300 = reshape(r_standard_old300, 1, n_sim);
r_comp_old300 = reshape(permute(r_comp_old300, [1 3 2]), cond, n_sim);


% young
for tt = 1:cond
    for ss = 1:n_sim
        for aff = 1:n_aff_old
            
            spk_o300{tt,ss,aff} = r_comp_old300{tt,ss}.responses(1,aff).spikes;
            spk_s_o300{ss,aff} =  r_standard_old300{ss}.responses(1,aff).spikes;
            
        end
    end
end


% reshape spike matrices
spike_o300 = reshape(permute(spk_o300,[2 1 3]), cond*n_sim, n_aff_old);

clear r_standard_old300 r_comp_old300

%% Elderly 30 Hz
load('FD_30Hz_old.mat');

% Extract spike trains
n_sim = 150;
cond = 5;

r_standard_old30 = reshape(r_standard_old30, 1, n_sim);
r_comp_old30 = reshape(permute(r_comp_old30, [1 3 2]), cond, n_sim);


% young
for tt = 1:cond
    for ss = 1:n_sim
        for aff = 1:n_aff_old
            
            spk_o30{tt,ss,aff} = r_comp_old30{tt,ss}.responses(1,aff).spikes;
            spk_s_o30{ss,aff} =  r_standard_old30{ss}.responses(1,aff).spikes;
            
        end
    end
end


% reshape spike matrices
spike_o30 = reshape(permute(spk_o30,[2 1 3]), cond*n_sim, n_aff_old);

save('spikes_old.mat','spk_s_o300','spike_o300','spk_s_o30','spike_o30')%'spk_s_y300','spike_y300','spk_s_y30','spike_y30',...
