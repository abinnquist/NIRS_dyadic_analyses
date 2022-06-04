%----------------
%PHASE SYNCHRONY 
%----------------
%instantaneous phase synchrony for all possible channel pairs - code from: Pedersen et al., 2017
currdir=cd;
codedir='C:\Users\Mike\Documents\MATLAB\';
%% Analysis to run
compile=0;

%% Choose either to compile data or load if already compiled
preprocess_dir=uigetdir('','Choose Data Directory');
addpath(genpath(strcat(currdir,filesep,'dependencies')))
addpath(genpath(strcat(codedir,'nirs-toolbox-master')))

dataprefix='SS';
ch_reject=3;
numScans=5;
zdim=1;
samprate = 3.9063; %change this to whatever sampling rate you had 
length_diss=1810; % Shortest conversation for both discussions

%% Load in the correct data 
if compile
    [deoxy3D,oxy3D]= compiledyadicNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans,zdim);

    %Save a .mat file to the preprocessing directory 
    save(strcat(preprocess_dir,filesep,'compiled_data'),'deoxy3D','oxy3D');
else
    load(strcat(preprocess_dir,filesep,'compiled_data'));
end

%% Run Phase synchrony 
%need a narrow frequency band for phase synchrony to work
[~, numchans, numdyads] = size(oxy3D(1).sub1);
dyads_phasesynch_fast1 = nan(numchans*2,numchans*2,numdyads); %slow-4 freq band
dyads_phasesynch_slow1 = nan(numchans*2,numchans*2,numdyads); %slow-5 freq band

dyads_phasesynch_fast2 = nan(numchans*2,numchans*2,numdyads); %slow-4 freq band
dyads_phasesynch_slow2 = nan(numchans*2,numchans*2,numdyads); %slow-5 freq band

num_good_channels = nan(1,numdyads);

for d=1:2 % number of convos
    for a=1:numdyads
        subj1 = oxy3D(d).sub1(1:length_diss,:,a);
        subj2 = oxy3D(d).sub2(1:length_diss,:,a);

        goodinds1 = find(~isnan(subj1(1,:))); 
        goodinds2 = find(~isnan(subj2(1,:)));
        goodindscheck1 = ~isnan(subj1(1,:));
        goodindscheck2 = ~isnan(subj2(1,:));
        goodindscheck = [goodindscheck1 goodindscheck2];
        subj1brief = subj1(:,goodinds1); %Remove NaN channels
        subj2brief = subj2(:,goodinds2);

        subj1filtered = nan(size(subj1brief));
        subj2filtered = nan(size(subj2brief));
        % Frequency to filter can be adjusted
        for x=1:size(subj1brief,2)
            subj1filtered(:,x) = hmrBandpassFilt(subj1brief(:,x), samprate, 0.03, 0.07);
        end
        for x=1:size(subj2brief,2)
            subj2filtered(:,x) = hmrBandpassFilt(subj2brief(:,x), samprate, 0.03, 0.07);
        end
        subjboth = [subj1filtered subj2filtered];
        [subjbothinnov,f] = nirs.math.innovations(subjboth,20); %removes auto-correlation
        num_good_channels(a)=size(subjboth,2);
        PSfast = angle(hilbert(subjbothinnov')); % instantaneous phases of channelxtime data
        PSfast_matrix = zeros(size(PSfast,1),size(PSfast,1),size(PSfast,2));
        for time_point = 1:size(PSfast, 2)
            %channel by channel phase synchrony at each time point - CxCxT
            PSfast_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSfast(:,time_point)', PSfast(:,time_point))));
        end
        for x=1:(size(subj1,2)*2)
            if ~goodindscheck(x)
                PSfast_matrix = [PSfast_matrix(:,1:(x-1),:) nan(size(PSfast_matrix,1),1,size(PSfast,2)) PSfast_matrix(:,x:end,:)];
            end
        end
        for x=1:(size(subj1,2)*2)
            if ~goodindscheck(x)
                PSfast_matrix = [PSfast_matrix(1:(x-1),:,:); nan(1,size(PSfast_matrix,2),size(PSfast,2)); PSfast_matrix(x:end,:,:)];
            end
        end
        
        if d==1
            dyads_phasesynch_fast1(:,:,a) = nanmean(PSfast_matrix,3);
        else
            dyads_phasesynch_fast2(:,:,a) = nanmean(PSfast_matrix,3);
        end

        subj1filtered = nan(size(subj1brief));
        subj2filtered = nan(size(subj2brief));
        for x=1:size(subj1brief,2)
            subj1filtered(:,x) = hmrBandpassFilt(subj1brief(:,x), samprate, 0.01, 0.03);
        end
        for x=1:size(subj2brief,2)
            subj2filtered(:,x) = hmrBandpassFilt(subj2brief(:,x), samprate, 0.01, 0.03);
        end
        subjboth = [subj1filtered subj2filtered];
        [subjbothinnov,f] = nirs.math.innovations(subjboth,20);
        num_good_channels(a)=size(subjboth,2);
        PSslow = angle(hilbert(subjbothinnov')); % instantaneous phases of channelxtime data
        PSslow_matrix = zeros(size(PSslow,1),size(PSslow,1),size(PSslow,2));
        for time_point = 1:size(PSslow, 2)
            %channel by channel phase synchrony at each time point - CxCxT
            PSslow_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSslow(:,time_point)', PSslow(:,time_point))));
        end
        for x=1:(size(subj1,2)*2)
            if ~goodindscheck(x)
                PSslow_matrix = [PSslow_matrix(:,1:(x-1),:) nan(size(PSslow_matrix,1),1,size(PSslow,2)) PSslow_matrix(:,x:end,:)];
            end
        end
        for x=1:(size(subj1,2)*2)
            if ~goodindscheck(x)
                PSslow_matrix = [PSslow_matrix(1:(x-1),:,:); nan(1,size(PSslow_matrix,2),size(PSslow,2)); PSslow_matrix(x:end,:,:)];
            end
        end

        if d==1
            dyads_phasesynch_slow1(:,:,a) = nanmean(PSslow_matrix,3);
        else
            dyads_phasesynch_slow2(:,:,a) = nanmean(PSslow_matrix,3);
        end
    end
end


%phase synch null dist: pairing random subjects who weren't in same dyad together.
iterations = 100; %number of bootstrapped null samples to take; small number here for speedier demo but should be 100+
newchans=numchans*2;
nulldist_fast = nan(newchans,newchans,iterations);
nulldist_slow = nan(newchans,newchans,iterations);

randsubj1 = randi([1 numdyads],1,iterations);
randsubj2 = randi([1 numdyads],1,iterations);
evenodd1 = randi([0 1],1,iterations);
evenodd2 = randi([0 1],1,iterations);

fprintf('\n\t Bootstrapping Phase Synchrony Null Distribution ...\n')
reverseStr = '';
Elapsedtime = tic;
for a=1:iterations
    if randsubj1(a)~=randsubj2(a)        
        if evenodd1(a)
            subj1 = oxy3D(d).sub1(1:length_diss,:,randsubj1(a));
        else
            subj1 = oxy3D(d).sub2(1:length_diss,:,randsubj1(a));
        end
        if evenodd2(a)
            subj2 = oxy3D(d).sub1(1:length_diss,:,randsubj2(a));
        else
            subj2 = oxy3D(d).sub2(1:length_diss,:,randsubj2(a));
        end
%         length1 = length(subj1);
%         length2 = length(subj2);
%         smallerlength = min(length1,length2);
%         subj1 = subj1(1:smallerlength,:);
%         subj2 = subj2(1:smallerlength,:);
    
        goodinds1 = find(~isnan(subj1(1,:)));
        goodinds2 = find(~isnan(subj2(1,:)));
        goodindscheck1 = ~isnan(subj1(1,:));
        goodindscheck2 = ~isnan(subj2(1,:));
        goodindscheck = [goodindscheck1 goodindscheck2];
        subj1brief = subj1(:,goodinds1);
        subj2brief = subj2(:,goodinds2);
    
        subj1fast = nan(size(subj1brief));
        subj2fast = nan(size(subj2brief));
        subj1slow = nan(size(subj1brief));
        subj2slow = nan(size(subj2brief));
    
        for x=1:size(subj1brief,2)
            subj1fast(:,x) = hmrBandpassFilt(subj1brief(:,x), samprate, 0.03, 0.07);
            subj1slow(:,x) = hmrBandpassFilt(subj1brief(:,x), samprate, 0.01, 0.027);
        end
        for x=1:size(subj2brief,2)
            subj2fast(:,x) = hmrBandpassFilt(subj2brief(:,x), samprate, 0.03, 0.07);
            subj2slow(:,x) = hmrBandpassFilt(subj2brief(:,x), samprate, 0.01, 0.027);
        end
    
        oxybothfast = [subj1fast subj2fast];
        oxybothslow = [subj1slow subj2slow];
    
        PSfast = angle(hilbert(oxybothfast')); % instantaneous phases of channelxtime data
        PSfast_matrix = zeros(size(PSfast,1),size(PSfast,1),size(PSfast,2));
        PSslow = angle(hilbert(oxybothslow')); 
        PSslow_matrix = zeros(size(PSslow,1),size(PSslow,1),size(PSslow,2));
    
        for time_point = 1:size(PSfast, 2)
        %channel by channel phase synchrony at each time point - CxCxT
            PSfast_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSfast(:,time_point)', PSfast(:,time_point))));
            PSslow_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSslow(:,time_point)', PSslow(:,time_point))));
        end
    
        for x=1:newchans
            if ~goodindscheck(x)
                PSfast_matrix = [PSfast_matrix(:,1:(x-1),:) nan(size(PSfast_matrix,1),1,size(PSfast,2)) PSfast_matrix(:,x:end,:)];
                PSslow_matrix = [PSslow_matrix(:,1:(x-1),:) nan(size(PSslow_matrix,1),1,size(PSslow,2)) PSslow_matrix(:,x:end,:)];
            end
        end
        for x=1:newchans
            if ~goodindscheck(x)
                PSfast_matrix = [PSfast_matrix(1:(x-1),:,:); nan(1,size(PSfast_matrix,2),size(PSfast,2)); PSfast_matrix(x:end,:,:)];
                PSslow_matrix = [PSslow_matrix(1:(x-1),:,:); nan(1,size(PSslow_matrix,2),size(PSslow,2)); PSslow_matrix(x:end,:,:)];
            end
        end
    
        nulldist_fast(:,:,a) = nanmean(PSfast_matrix,3);
        nulldist_slow(:,:,a) = nanmean(PSslow_matrix,3);
    else
        nulldist_fast(:,:,a) = nan(newchans,newchans);
        nulldist_slow(:,:,a) = nan(newchans,newchans);
    end
    msg = sprintf('\n\t iteration %d/%d ...',a,iterations);
    fprintf([reverseStr,msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end
Elapsedtime = toc(Elapsedtime);
fprintf('\n\t Elapsed time: %g seconds ...\n', Elapsedtime);

nulldist_fast_inter = nulldist_fast(1:numchans,numchans+1:end,:);
nulldist_slow_inter = nulldist_slow(1:numchans,numchans+1:end,:);
nulldist_fast_inter = nulldist_fast_inter(~isnan(nulldist_fast_inter));%makes a single vector of inter connections, ignores channel label
nulldist_slow_inter = nulldist_slow_inter(~isnan(nulldist_slow_inter));
fisherz = @(r)(log(1+r)-log(1-r))/2;
nulldist_fast_z = arrayfun(fisherz, nulldist_fast_inter);
nulldist_slow_z = arrayfun(fisherz, nulldist_slow_inter);
slst_fast = sort(nulldist_fast_z);
slst_slow = sort(nulldist_slow_z);
slength_fast = length(slst_fast);
slength_slow = length(slst_slow);
phasecutoff_fast = slst_fast(end-round(slength_fast*5/100));
phasecutoff_slow = slst_slow(end-round(slength_slow*5/100));
    %For both Convos
dyad_fast_z1 = arrayfun(fisherz, dyads_phasesynch_fast1(1:numchans,numchans+1:end,:));
dyad_slow_z1 = arrayfun(fisherz, dyads_phasesynch_slow1(1:numchans,numchans+1:end,:));
dyad_fast_z2 = arrayfun(fisherz, dyads_phasesynch_fast2(1:numchans,numchans+1:end,:));
dyad_slow_z2 = arrayfun(fisherz, dyads_phasesynch_slow2(1:numchans,numchans+1:end,:));

%network estimation based on similarity matrices and null dist made above
binarynetwork_fast1 = zeros(numchans,numchans,numdyads);
binarynetwork_slow1 = zeros(numchans,numchans,numdyads);
binarynetwork_fast2 = zeros(numchans,numchans,numdyads);
binarynetwork_slow2 = zeros(numchans,numchans,numdyads);
for z=1:numdyads
binarynetwork_fast1(:,:,z) = dyad_fast_z1(:,:,z)>phasecutoff_fast;
binarynetwork_slow1(:,:,z) = dyad_slow_z1(:,:,z)>phasecutoff_slow;

binarynetwork_fast2(:,:,z) = dyad_fast_z2(:,:,z)>phasecutoff_fast;
binarynetwork_slow2(:,:,z) = dyad_slow_z2(:,:,z)>phasecutoff_slow;
end

%Find p-values of all channel pairings compared to the null dist
pvals_fast1 = ones(numchans,numchans);
pvals_slow1 = ones(numchans,numchans);
pvals_fast2 = ones(numchans,numchans);
pvals_slow2 = ones(numchans,numchans);
for x=1:numchans
    for y=1:numchans
        fastval1 = nanmean(dyad_fast_z1(x,y,:));
        slowval1 = nanmean(dyad_slow_z1(x,y,:));
        fastval2 = nanmean(dyad_fast_z2(x,y,:));
        slowval2 = nanmean(dyad_slow_z2(x,y,:));
        
        fast_pval1 = slst_fast>fastval1;
        slow_pval1 = slst_slow>slowval1;
        fast_pval2 = slst_fast>fastval2;
        slow_pval2 = slst_slow>slowval2;
        
        pvals_fast1(x,y)=sum(fast_pval1)/slength_fast;
        pvals_slow1(x,y)=sum(slow_pval1)/slength_slow;
        pvals_fast2(x,y)=sum(fast_pval2)/slength_fast;
        pvals_slow2(x,y)=sum(slow_pval2)/slength_slow;
    end
end