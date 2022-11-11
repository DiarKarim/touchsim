
%% Frequency Discrimination - Simulate response
% Young vs Elderly - 2 point discrimination
% Young's Modulus --> for young adult is 50kPa
%                 --> for elderly is 35kPa

clear 
addpath /Users/wingam/Downloads/touchSim-master
cd /Users/wingam/Downloads/touchSim-master

% create afferent population
% young

[xRA,yRA] = meshgrid(-5:0.8:6,-5:0.8:6); % density for young
[xSA1,ySA1] = meshgrid(-5:1.1:6,-5:1.1:6);
[xPC,yPC] = meshgrid(-5:2:6,-5:2:6);

sdRA = 0.8 * 0.2;
sdSA1 = 1.1 * 0.2;
sdPC = 2 * 0.2;
% 
for tt = 1:length(xRA)
    for ii = 1:length(xRA)
        xRA(tt,ii)  = xRA(tt,ii) + (0 + sdRA.*randn(1,1));
        yRA(tt,ii)  = yRA(tt,ii) + (0 + sdRA.*randn(1,1));

    end
end

for tt = 1:length(xSA1)
    for ii = 1:length(xSA1)
        xSA1(tt,ii)  = xSA1(tt,ii) + (0 + sdSA1.*randn(1,1));
        ySA1(tt,ii)  = ySA1(tt,ii) + (0 + sdSA1.*randn(1,1));

    end
end

for tt = 1:length(xPC)
    for ii = 1:length(xPC)
        xPC(tt,ii)  = xPC(tt,ii) + (0 + sdPC.*randn(1,1));
        yPC(tt,ii)  = yPC(tt,ii) + (0 + sdPC.*randn(1,1));

    end
end

aff_y = AfferentPopulation();
aff_y.add_afferents('SA1',[xSA1(:) ySA1(:)]);
aff_y.add_afferents('RA',[xRA(:) yRA(:)]);
aff_y.add_afferents('PC',[xPC(:) yPC(:)]);

n_SA1_y = sum(aff_y.iSA1);
n_RA_y = sum(aff_y.iRA);
n_PC_y = sum(aff_y.iPC);
n_aff_y = n_SA1_y + n_RA_y + n_PC_y;

figure(1)
plot(aff_y,[],'region','D2d')

[xRA,yRA] = meshgrid(-5:1.4:6,-5:1.4:6); % density for elderly
[xSA1,ySA1] = meshgrid(-5:1.8:6,-5:1.8:6);
[xPC,yPC] = meshgrid(-5:2:6,-5:2:6);

sdRA = 1.4 * 0.2;
sdSA1 = 1.8 * 0.2;
sdPC = 2 * 0.2;

for tt = 1:length(xRA)
    for ii = 1:length(xRA)
        xRA(tt,ii)  = xRA(tt,ii) + (0 + sdRA.*randn(1,1));
        yRA(tt,ii)  = yRA(tt,ii) + (0 + sdRA.*randn(1,1));

    end
end

for tt = 1:length(xSA1)
    for ii = 1:length(xSA1)
        xSA1(tt,ii)  = xSA1(tt,ii) + (0 + sdSA1.*randn(1,1));
        ySA1(tt,ii)  = ySA1(tt,ii) + (0 + sdSA1.*randn(1,1));

    end
end

for tt = 1:length(xPC)
    for ii = 1:length(xPC)
        xPC(tt,ii)  = xPC(tt,ii) + (0 + sdPC.*randn(1,1));
        yPC(tt,ii)  = yPC(tt,ii) + (0 + sdPC.*randn(1,1));

    end
end

aff_old = AfferentPopulation();
aff_old.add_afferents('SA1',[xSA1(:) ySA1(:)]);
aff_old.add_afferents('RA',[xRA(:) yRA(:)]);
aff_old.add_afferents('PC',[xPC(:) yPC(:)]);

n_SA1_old = sum(aff_old.iSA1);
n_RA_old = sum(aff_old.iRA);
n_PC_old = sum(aff_old.iPC);
n_aff_old = n_SA1_old + n_RA_old + n_PC_old;
% 
% figure(2)
%plot(aff_old,[],'region','D2d')
% %save('afferent','aff_y','aff_old')

%load('afferent','aff_y','aff_old')

%% Create standard stimulus - 30 Hz
% ---------------------------------------------------------------------
% ------------------------------ Standard stimulus --------------------
% ---------------------------------------------------------------------

% This is for the actuator 
sf = 5000; % Sampling frequency (samples per second)
len = 0.5; % seconds
freq = 30; %Hz
rad = 5; % mm
amp = [0.4 0.6 0.8]; % mm 3 amps for varibility so that the model responds more realistically 
ramp_len = 0.1; % Ramp up and down i.e. 100 ms regardless of the stimulus duration 

xy = [0 0]; % pin coordinate i.e. vibrotactor location at the center of the finger 

for ii = 1:length(amp)
stim = stim_sine(freq,amp(ii),0,len,xy,sf,ramp_len,rad); % This creates the signal 

trace = stim.trace; % i.e. trace of the sine wave 

standardStim{ii} = Stimulus(trace,xy,5000,rad); % generate stimulus object i.e. apply the signal
% 
% figure(ii)
% plot(standardStim{ii})
end
% create comparison stimuli

freq = 30;
deltaFreq = freq + [1 2 3 5 7];

cond = length(deltaFreq);

for tt = 1:cond
    for ii = 1:length(amp)
        
        c_stim = stim_sine(deltaFreq(tt),amp(ii),0,len,xy,sf,ramp_len,rad);
        
        trace = c_stim.trace;
        
        comparisonStim{tt,ii} = Stimulus(trace,xy,5000,rad); % generate stimulus object
        
        %plot(comparisonStim{4,2}.trace)
        
    end
end

% Stimulus function also contains a definition for skin properties (youngs
% modulus and poisson ration). At the moment we only changed the youngs
% modulus but this has to be done in the source code ()

% Simulate response
n_sim = 50; % Number of trials to generate 
r_standard_y30 = cell(length(amp), n_sim); % Declare variable 

% Diar -  Generate response
% for standard frequency
for ii = 1:length(amp)  
    for tt = 1:n_sim
        % young
        r_standard_y30{ii,tt} = aff_y.response(standardStim{ii});
        
        %     rates_s_y{tt} = r_single_y{tt}.rate;
        %     rates_s_y{tt}(a.iSA1)=rates_s_y{tt}(a.iSA1)/max(rates_s_y{tt}(a.iSA1));
        %     rates_s_y{tt}(a.iRA)=rates_s_y{tt}(a.iRA)/max(rates_s_y{tt}(a.iRA));
        %     rates_s_y{tt}(a.iPC)=rates_s_y{tt}(a.iPC)/max(rates_s_y{tt}(a.iPC));
    end
end

% ---------------------------------------------------------------------
% ------------------------------ 30 Hz stimulus -----------------------
% ---------------------------------------------------------------------
%
figure(1)
subplot(1,2,1)
plot(r_standard_y30{1})
subplot(1,2,2)
plot(r_standard_y30{2})

%
r_comp_y30 = cell(length(deltaFreq),length(amp), n_sim);
 
% for comparison frequency
for jj = 1:length(deltaFreq)
    for ii = 1:length(amp)
        for tt = 1:n_sim
            
            r_comp_y30{jj,ii,tt} = aff_y.response(comparisonStim{jj,ii});
        end
    end
end



%
% figure(1)
% plot(r_single_y{1,3})
% figure(2)
% plot(r_young{6,1})
% plot(a,[],'rate',rates_y{4},'region','D2')
% figure(6)
% plot(b,[],'rate',rates_old{4},'region','D2')
% figure(7)
% plot(a,[],'rate',rates_s_y{4},'region','D2')
% figure(8)
% plot(b,[],'rate',rates_s_o{4},'region','D2')

cd /Users/d.deflorio/Desktop/simulation/FrequencyDiscrimination

save('FD_30Hz_young.mat','r_comp_y30','r_standard_y30','-v7.3')

%% Create standard stimulus - 300 Hz
% ---------------------------------------------------------------------
% ------------------------------ 300 Hz stimulus -----------------------
% ---------------------------------------------------------------------

sf = 10000; % Sampling frequency (samples per second)
len = 0.5; % seconds
freq = 300; %Hz
rad = 5; % mm
amp = [0.4 0.6 0.8]; % mm
ramp_len = 0.1;

xy = [0 0]; % pin coordinate

for ii = 1:length(amp)
stim = stim_sine(freq, amp(ii), 0, len, xy, sf, ramp_len, rad);

trace = stim.trace;

standardStim{ii} = Stimulus(trace,xy,sf,rad); % generate stimulus object
% 
% figure(ii)
% plot(standardStim{ii})
end
% create comparison stimuli

freq = 300;
deltaFreq = freq + [2 5 10 18 30];

cond = length(deltaFreq);

for tt = 1:cond
    for ii = 1:length(amp)
        
        c_stim = stim_sine(deltaFreq(tt),amp(ii),0,len,xy,sf,ramp_len,rad);
        
        trace = c_stim.trace;
        
        comparisonStim{tt,ii} = Stimulus(trace,xy,sf,rad); % generate stimulus object
        
        %plot(comparisonStim{4,2}.trace)
        
    end
end

% Simulate response
n_sim = 50;

r_standard_y300 = cell(length(amp), n_sim);

% for standard frequency
for ii = 1:length(amp)  
    for tt = 1:n_sim
        % young
        r_standard_y300{ii,tt} = aff_y.response(standardStim{ii});
        
        %     rates_s_y{tt} = r_single_y{tt}.rate;
        %     rates_s_y{tt}(a.iSA1)=rates_s_y{tt}(a.iSA1)/max(rates_s_y{tt}(a.iSA1));
        %     rates_s_y{tt}(a.iRA)=rates_s_y{tt}(a.iRA)/max(rates_s_y{tt}(a.iRA));
        %     rates_s_y{tt}(a.iPC)=rates_s_y{tt}(a.iPC)/max(rates_s_y{tt}(a.iPC));
    end
end

%
% figure(1)
% subplot(1,2,1)
% plot(r_standard_y300{1})
% subplot(1,2,2)
% plot(r_standard_y300{2})

r_comp_y300 = cell(length(deltaFreq),length(amp), n_sim);

% for comparison frequency
for jj = 1:length(deltaFreq)
    for ii = 1:length(amp)
        for tt = 1:n_sim
            
            r_comp_y300{jj,ii,tt} = aff_y.response(comparisonStim{jj,ii});
        end
    end
end



%
% figure(1)
% plot(r_single_y{1,3})
% figure(2)
% plot(r_young{6,1})
% plot(a,[],'rate',rates_y{4},'region','D2')
% figure(6)
% plot(b,[],'rate',rates_old{4},'region','D2')
% figure(7)
% plot(a,[],'rate',rates_s_y{4},'region','D2')
% figure(8)
% plot(b,[],'rate',rates_s_o{4},'region','D2')

%cd /Users/d.deflorio/Desktop/simulation/FrequencyDiscrimination

save('FD_300Hz_young.mat','r_comp_y300','r_standard_y300','-v7.3')


%% FOR ELDERLY - low aff, high stiff


%% Create standard stimulus - 30 Hz

sf = 5000; % Sampling frequency (samples per second)
len = 0.5; % seconds
freq = 30; %Hz
rad = 5; % mm
amp = [0.4 0.6 0.8]; % mm
ramp_len = 0.1;

xy = [0 0]; % pin coordinate

for ii = 1:length(amp)
stim = stim_sine(freq,amp(ii),[],len,xy,sf,ramp_len,rad);

trace = stim.trace;

standardStim{ii} = Stimulus_Old(trace,xy,5000,rad); % generate stimulus object
% 
% figure(ii)
% plot(standardStim{ii})
end
%% create comparison stimuli

freq = 30;
deltaFreq = freq + [1 2 3 5 7];

cond = length(deltaFreq);

for tt = 1:cond
    for ii = 1:length(amp)
        
        c_stim = stim_sine(deltaFreq(tt),amp(ii),[],len,xy,sf,ramp_len,rad);
        
        trace = c_stim.trace;
        
        comparisonStim{tt,ii} = Stimulus_Old(trace,xy,5000,rad); % generate stimulus object
        
        %plot(comparisonStim{4,2}.trace)
        
    end
end

%% Simulate response
n_sim = 50;

r_standard_old30 = cell(length(amp), n_sim);

% for standard frequency
for ii = 1:length(amp) 
    for tt = 1:n_sim
        % young
        r_standard_old30{ii,tt} = aff_old.response(standardStim{ii});
        
        %     rates_s_y{tt} = r_single_y{tt}.rate;
        %     rates_s_y{tt}(a.iSA1)=rates_s_y{tt}(a.iSA1)/max(rates_s_y{tt}(a.iSA1));
        %     rates_s_y{tt}(a.iRA)=rates_s_y{tt}(a.iRA)/max(rates_s_y{tt}(a.iRA));
        %     rates_s_y{tt}(a.iPC)=rates_s_y{tt}(a.iPC)/max(rates_s_y{tt}(a.iPC));
    end
end

%%
figure(1)
subplot(1,2,1)
plot(r_standard_y30{1})
subplot(1,2,2)
plot(r_standard_old30{2})

%%
r_comp_old30 = cell(length(deltaFreq),length(amp), n_sim);

% for comparison frequency
for jj = 1:length(deltaFreq)
    for ii = 1:length(amp)
        for tt = 1:n_sim
            
            r_comp_old30{jj,ii,tt} = aff_old.response(comparisonStim{jj,ii});
        end
    end
end



%
% figure(1)
% plot(r_single_y{1,3})
% figure(2)
% plot(r_young{6,1})
% plot(a,[],'rate',rates_y{4},'region','D2')
% figure(6)
% plot(b,[],'rate',rates_old{4},'region','D2')
% figure(7)
% plot(a,[],'rate',rates_s_y{4},'region','D2')
% figure(8)
% plot(b,[],'rate',rates_s_o{4},'region','D2')

%cd /Users/d.deflorio/Desktop/simulation/FrequencyDiscrimination

save('FD_30Hz_old.mat','r_comp_old30','r_standard_old30','-v7.3')

%% Create standard stimulus - 300 Hz

sf = 10000; % Sampling frequency (samples per second)
len = 0.5; % seconds
freq = 300; %Hz
rad = 5; % mm
amp = [0.4 0.6 0.8]; % mm
ramp_len = 0.1;

xy = [0 0]; % pin coordinate

for ii = 1:length(amp)
stim = stim_sine(freq,amp(ii),[],len,xy,sf,ramp_len,rad);

trace = stim.trace;

standardStim{ii} = Stimulus_Old(trace,xy,sf,rad); % generate stimulus object
% 
% figure(ii)
% plot(standardStim{ii})
end
%% create comparison stimuli

freq = 300;
deltaFreq = freq + [2 5 10 18 30];

cond = length(deltaFreq);

for tt = 1:cond
    for ii = 1:length(amp)
        
        c_stim = stim_sine(deltaFreq(tt),amp(ii),[],len,xy,sf,ramp_len,rad);
        
        trace = c_stim.trace;
        
        comparisonStim{tt,ii} = Stimulus_Old(trace,xy,sf,rad); % generate stimulus object
        
        %plot(comparisonStim{4,2}.trace)
        
    end
end

%% Simulate response
n_sim = 50;
r_standard_old300 = cell(length(amp), n_sim);

% for standard frequency
for ii = 1:length(amp)   
    for tt = 1:n_sim
        % old
        r_standard_old300{ii,tt} = aff_old.response(standardStim{ii});
        
        %     rates_s_o{tt} = r_single_o{tt}.rate;
        %     rates_s_o{tt}(a.iSA1)=rates_s_o{tt}(a.iSA1)/max(rates_s_o{tt}(a.iSA1));
        %     rates_s_o{tt}(a.iRA)=rates_s_o{tt}(a.iRA)/max(rates_s_o{tt}(a.iRA));
        %     rates_s_o{tt}(a.iPC)=rates_s_o{tt}(a.iPC)/max(rates_s_o{tt}(a.iPC));
    end
end

%%
% figure(1)
% subplot(1,2,1)
% plot(r_standard_o300{1})
% subplot(1,2,2)
% plot(r_standard_o300{2})

%%
r_comp_old300 = cell(length(deltaFreq), length(amp), n_sim);

% for comparison frequency
for jj = 1:length(deltaFreq)
    for ii = 1:length(amp)
        for tt = 1:n_sim
            
            r_comp_old300{jj,ii,tt} = aff_old.response(comparisonStim{jj,ii});
        end
    end
end

%
% figure(1)
% plot(r_single_y{1,3})
% figure(2)
% plot(r_young{6,1})
% plot(a,[],'rate',rates_y{4},'region','D2')
% figure(6)
% plot(b,[],'rate',rates_old{4},'region','D2')
% figure(7)
% plot(a,[],'rate',rates_s_y{4},'region','D2')
% figure(8)
% plot(b,[],'rate',rates_s_o{4},'region','D2')

%cd /Users/d.deflorio/Desktop/simulation/FrequencyDiscrimination

save('FD_300Hz_old.mat','r_comp_old300','r_standard_old300','-v7.3')
