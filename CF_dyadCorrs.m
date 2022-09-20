clc; clear
%% Instructions %%
% Run from the main folder 'NIRS_dyadic_analysis'
% Dependencies for analyses located in helperFuncs
% Dependencies for imaging: spm12, xjview, and mni_coordinates (in folder)
%
% Data must be preprocessed before using this script (see
% https://github.com/abinnquist/fNIRSpreProcessing)
%
% Before running the script make sure to choose the analyses you want to run 
% by turning on (1) or turning off (0) the tasks you want to be performed. 
% You can change the properties for a more strict/leneint FDR cutoff or
% different trim time.
%
% If you prefer to specify your data path instead of selecting each time
% comment out the uigetdir command and uncomment the command above it. Make
% sure to specify the location of the preprocessed NIRS data
%% Analyses to run
% Set to zero if you do not want to perform the task
compile=0;      % If the data has yet to be compiled after preprocessing 
oxyOnly=1;      % 0=z_deoxy; 1=z_oxy;
chCorr=0;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=1;     % Dyadic correlation over entire conversation per area of the brain
FDR=1;          % False discovery rate correction
writeXL=1;      % Write the data to an excel sheet(s) in the preprocess_dir
image=0;        % Will prompt you in command window for mni, conversation & dyad 

%% Load directory & varibles to use
% You can set the directory w/ line 3 or use the get directory in line 4
%preprocess_dir= ''; %Specify path instead of selecting each run
preprocess_dir=uigetdir('','Choose Data Directory');
addpath(genpath('helperFuncs'))
addpath(genpath('dependencies'))

%% Properties
% These can be changed based on the study and what FDR cutoff is prefered
dataprefix='0';
ch_reject=3; % 1=no channel rejection, 2=reject noisy, 3=reject noisy & uncertain
cutoff=0.1; % Cut-off p-value to use for FDR (false discovery rate) 
length_scan=2442; % Frame to trim all subjects (shortest convo)
numdyads=54; % Number of dyads
numchans=42; % Number of channels
numareas=4;  % 4=all ch's mPFC, lPFC, TPJ, & VIS; 5=subsets

%% Compile data for easier analysis
% This will compile the two conversations of interest. Creates twelve 3D
% matrices (time,channels,dyad) based on type of oxy (deoxy & oxy),  
% discussion (affiliation and conflict), and subject (a and b), (2x2x2=8).
if compile
    numScans=2;
    [deoxy3D,oxy3D]= compiledyadicNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans)

    %Save a .mat file to the preprocessing directory 
    save(strcat(preprocess_dir,filesep,'compiled_data'),'deoxy3D','oxy3D');

    clearvars -except preprocess_dir numdyads numchans numareas length_scan oxyOnly chCorr areaCorr cutoff FDR writeXL image
end

%% Mean Timecourse Synchrony per ROI per Dyad %%
if areaCorr
    montageMatch=readtable(strcat(preprocess_dir,filesep,'montageMatch.csv'));

    if numareas==4 
        areas={7:14;[1:6,15:20];[25:30,36:40];[21:24,32:35]}; %mPFC, lPFC, TPJ, VM
        areas(:,2)={[7,11:14];[3:6,15:16,19:20];[25:26,28,30,36:37,39,41];[21,23:24,32,34:35]}; %mPFC, lPFC, TPJ, VM
    elseif numareas==5
        areas={7:14;[1:6,15:20];[25:30,36:40];[21:24,32:35];[31,42]}; %mPFC, lPFC, TPJ, VM, FF
        areas(:,2)={[7,11:14];[3:6,15:16,19:20];[25:26,28,30,36:37,39,41];[21,23:24,32,34:35];[31,42]}; %mPFC, lPFC, TPJ, VM, FF
    end  
   
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'oxy3D')
        for dy=1:numdyads
            row1=montageMatch.Montage1(dy);
            row2=montageMatch.Montage2(dy);
    
            for ar=1:numareas
                z1_aff_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_scan,cell2mat(areas(ar,row1)),dy),2);
                z2_aff_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_scan,cell2mat(areas(ar,row2)),dy),2);
                z1_con_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_scan,cell2mat(areas(ar,row1)),dy),2);
                z2_con_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_scan,cell2mat(areas(ar,row2)),dy),2);
        
                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas(ar,row2)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas(ar,row2)),dy))));
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'deoxy3D')
        
        for dy=1:numdyads
            row1=montageMatch.Montage1(dy);
            row2=montageMatch.Montage2(dy);
    
            for ar=1:numareas
                z1_aff_areas(:,ar,dy) = nanmean(deoxy3D(1).sub1(1:length_scan,cell2mat(areas(ar,row1)),dy),2);
                z2_aff_areas(:,ar,dy) = nanmean(deoxy3D(1).sub2(1:length_scan,cell2mat(areas(ar,row2)),dy),2);
                z1_con_areas(:,ar,dy) = nanmean(deoxy3D(2).sub1(1:length_scan,cell2mat(areas(ar,row1)),dy),2);
                z2_con_areas(:,ar,dy) = nanmean(deoxy3D(2).sub2(1:length_scan,cell2mat(areas(ar,row2)),dy),2);
        
                missN_s1(dy,ar,1)=numel(find(isnan(deoxy3D(1).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(deoxy3D(1).sub2(1,cell2mat(areas(ar,row2)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(deoxy3D(2).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(deoxy3D(2).sub2(1,cell2mat(areas(ar,row2)),dy))));
            end
        end
    end

    %Creates a mask for missing channels for both subjects within each dyad
    loss=0.25;
    nmask_s1=zeros(numdyads,numareas,2);
    nmask_s2=zeros(numdyads,numareas,2);
    for ds=1:2
        for dy=1:numdyads
            row1=montageMatch.Montage1(dy);
            row2=montageMatch.Montage2(dy);
            for ar=1:numareas
                if missN_s1(dy,ar,ds)>=loss*length(cell2mat(areas(ar,row1)))
                    nmask_s1(dy,ar,ds)=1;
                end
                if missN_s2(dy,ar,ds)>=loss*length(cell2mat(areas(ar,row2)))
                    nmask_s2(dy,ar,ds)=1;
                end
            end
        end
    end   
    nmask=nmask_s1+nmask_s2;
    nmask(nmask(:,:,:)==1)=0;
    nmask(nmask(:,:,:)==2)=1;
    nmask=logical(nmask);

    % Correlation of ROI x dyad x convo
    for area=1:numareas
        for dyad=1:numdyads
            a = z1_aff_areas(:,area,dyad);
            b = z2_aff_areas(:,area,dyad);
            [r_aff_areas(dyad,area),p_aff_areas(dyad,area)] = corr(a,b);

            c = z1_con_areas(:,area,dyad);
            d = z2_con_areas(:,area,dyad);
            [r_con_areas(dyad,area),p_con_areas(dyad,area)] = corr(c,d);
        end    
    end

    [mask_aff, ~]=fdr_bky(p_aff_areas,cutoff,'yes'); %Creates a mask for sig. areas
    [mask_con, ~]=fdr_bky(p_con_areas,cutoff,'yes'); %Creates a mask for sig. areas

    Sig_r_affiliation=r_aff_areas;
    Sig_r_affiliation(~mask_aff)=NaN;
    Sig_r_aff=Sig_r_affiliation;
    Sig_r_aff(~nmask(:,:,1))=NaN;
    Sig_r_affiliation(:,:,2)=Sig_r_aff;
    Sig_r_affiliation(:,:,3)=r_aff_areas;

    Sig_r_conflict=r_con_areas;
    Sig_r_conflict(~mask_con)=NaN;
    Sig_r_con=Sig_r_conflict;
    Sig_r_con(nmask(:,:,2))=NaN;
    Sig_r_conflict(:,:,2)=Sig_r_con;
    Sig_r_conflict(:,:,3)=r_con_areas;
     
    if writeXL
        %Sig_r_conflict=atanh(Sig_r_conflict);%z-scored r-values
        load(strcat('CF_NIRS',filesep,'CF_dyads.mat'))
        
        if oxyOnly
            typeOxy='oxy_'; %file naming for type of oxy
        else
            typeOxy='deoxy_';
        end
            
        if numareas==4 %variable names and file name for # of areas
            areas=["Subs","Dyads","mPFC","lPFC","TPJ","VIS"];
            aName='Four_';
        elseif numareas==5
            areas=["Subs","Dyads","mPFC","lPFC","TPJ","VMC","VIS"];
            aName='new_';
        end
        
        Sig_r_affiliation = repelem(Sig_r_affiliation,2,1);
        Sig_r_conflict = repelem(Sig_r_conflict,2,1);
        dyads=repelem(dyads,2,1);
        subs=repmat(1:2,1,numdyads)';

        Sig_r_aff1=array2table([subs,dyads,Sig_r_affiliation(:,:,1)],'VariableNames',areas);
        Sig_r_aff2=array2table([subs,dyads,Sig_r_affiliation(:,:,2)],'VariableNames',areas);
        r_values_aff=array2table([subs,dyads,Sig_r_affiliation(:,:,3)],'VariableNames',areas);

        Sig_r_con1=array2table([subs,dyads,Sig_r_conflict(:,:,1)],'VariableNames',areas);
        Sig_r_con2=array2table([subs,dyads,Sig_r_conflict(:,:,2)],'VariableNames',areas);
        r_values_con=array2table([subs,dyads,Sig_r_conflict(:,:,3)],'VariableNames',areas);
        
        %File name based on FDR used
        xlName=strcat(preprocess_dir,filesep,aName,typeOxy,'con_mask_',num2str(cutoff),'.xls'); 
        writetable(r_values_con,xlName,'sheet','r_vals')  %r-vals no mask or FDR
        writetable(Sig_r_con1,xlName,'sheet','r_vals_FDR')  %r-vals lost chans masked
        writetable(Sig_r_con2,xlName,'sheet','r_vals_FDRchLoss') %r-vals FDR & mask
    end
    save(strcat(preprocess_dir,filesep,aName,typeOxy,'CF_areas'),'Sig_r_conflict','Sig_r_affiliation')
    clearvars -except preprocess_dir numchans numdyads image
end

%% Image NIRS results
% The mni coordinates should be saved in the study folder or you can select
% wherever you have it saved
if image
    mniQ = 'Do you have the mni coordinates? 0=No, 1=Yes. \n';
    mni = input(mniQ);
    if mni==1
        [mniCds, mniPath] = uigetfile('*.mat','Choose the MNI coordinates .mat');
        mniCoords = strcat(mniPath,mniCds);
        mniCoords=struct2array(load(mniCoords));
        
        fileMade=imageNIRSvals(mniCoords); % Helper function in cd
        disp(fileMade)
    else
        disp('Must have an mni.mat to image data');
    end
end
clear
