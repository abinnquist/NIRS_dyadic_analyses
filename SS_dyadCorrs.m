clc;  clear
%% Instructions %%
% Run from the main folder 'NIRS_dyadic_analysis'
% Data must be preprocessed before using this script (see
% https://github.com/abinnquist/fNIRSpreProcessing)
%
% Dependencies for analyses located in helperFuncs
%
% Before running the script make sure to choose the analyses you want to run 
% by turning on (1), or turning off (0), the tasks you want to perform. 
% You can change the properties for a different trim time, more 
% strict/leneint FDR cutoff and/or channel rejection based on the three 
% outputs from preprocessing.
%
% If you prefer to specify your data path instead of selecting each time
% comment out the uigetdir command and uncomment the command above it. Make
% sure to specify the location of the preprocessed NIRS data
%% Analyses to run
% Set to zero if you do not want to perform the task
compile=0;      % If the data has yet to be compiled into all dyad matices
oxyOnly=1;      % 0=z_deoxy; 1=z_oxy
chCorr=0;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=1;     % Dyadic correlation over entire conversation per area of the brain
FDR=1;          % False discovery rate correction
writeXL=1;      % If you want to write the data to an excel sheet(s)

%% Set the Directory
% You can set the directory w/ line 3 or use the get directory in line 4
% preprocess_dir= ''; %Make sure to specify correct location
preprocess_dir=uigetdir('','Choose Data Directory');
addpath(genpath('helperFuncs'))

%% Properties
% These can be changed based on the study and what FDR cutoff is prefered
dataprefix='SS';
ch_reject=3;  % 1=no channel rejection, 2=reject noisy, 3=reject noisy & uncertain
cutoff=0.1; % Cut-off p-value to use for FDR (false discovery rate) 
length_diss1=1853; % Can be changed to the shortest conversation
length_diss2=1810; % Can be changed to the shortest conversation
numdyads=52; % Number of dyads
numchans=42; % Number of channels
numareas=5; % 4=mPFC, lPFC, TPJ, VIS; 5=mPFC, lPFC, SMS, TPJ & VIS; 6=vmPFC, dmPFC, vlPFC, SMS, TPJ
zdim=1; % 0=compile non-zscored; 1=compile z-scored

%% Compile subject data %%
% Compiles the two conversations of interest. Creates a 3D
% matrix for type of oxy (deoxy & oxy), discussion (1 and 2),  
% and subject (1 and 2) for eight (2x2x2=12) 3D matrices (time,channels,dyad)
if compile
    numScans=5;
    [deoxy3D,oxy3D]= compiledyadicNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans,zdim);

    %Save a .mat file to the preprocessing directory 
    save(strcat(preprocess_dir,filesep,'compiled_data'),'deoxy3D','oxy3D');
    
    clearvars -except preprocess_dir numdyads numchans numareas length_diss1 length_diss2 oxyOnly chCorr areaCorr cutoff FDR writeXL
end

%% Mean Timecourse Synchrony per ROI per Dyad %%
% Average activiation in area specific channels, then computes dyadic 
% correlation of those areas of the brain 
if areaCorr    
    montageMatch=readtable(strcat(preprocess_dir,filesep,'montageMatch.csv'));

    if numareas==4 
        areas={[1:3,5,8,10:12];[4,6:7,9,30:31];[17:21,23,36:40,42];[22,41]}; %mPFC, lPFC, TPJ, VIS
        areas(:,2)={[1:3,5,8,10];[4,6:7,9,30:31];[16:21,35:40];41:42}; %mPFC, lPFC, TPJ, VIS
        areas(:,3)={[1:3,5,8];[4,7,30:31];[14:18,33:37];[21:23,40:42]}; %mPFC, lPFC, TPJ, VIS
    elseif numareas==5
        areas={[1:3,5,8,10:12];[4,6:7,9,30:31];[14:16,26:29,33:35];[17:21,23,36:40,42];[22,41]}; %mPFC, lPFC, SSC, TPJ, VIS
        areas(:,2)={[1:3,5,8,10];[4,6:7,9,30:31];[13:14,24:25,32:34];[16:21,27,29,35:40];[22:23,41:42]}; %mPFC, lPFC, SSC, TPJ, VIS
        areas(:,3)={[1:3,5,8];[4,7,30:31];[13,24:25,32];[14:18,33:37];[21:23,40:42]}; %mPFC, lPFC, SSC, TPJ, VIS
    end
    
    z1_diss1_areas=nan(length_diss1,numareas,numdyads);
    z2_diss1_areas=nan(length_diss1,numareas,numdyads);
    z1_diss2_areas=nan(length_diss2,numareas,numdyads);
    z2_diss2_areas=nan(length_diss2,numareas,numdyads);
    
    missN_s1=zeros(numdyads,numareas,2);
    missN_s2=zeros(numdyads,numareas,2);
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'compiled_data'),'oxy3D');
        for dy=1:numdyads
            row1=montageMatch.Montage1(dy);
            row2=montageMatch.Montage2(dy);
    
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas(ar,row1)),dy),2);
                z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas(ar,row2)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas(ar,row1)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas(ar,row2)),dy),2);
        
                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas(ar,row2)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas(ar,row2)),dy))));
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'deoxy3D')
        for dy=1:numdyads
            row1=montageMatch.Montage1(dy);
            row2=montageMatch.Montage2(dy);
    
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(deoxy3D(1).sub1(1:length_diss1,cell2mat(areas(ar,row1)),dy),2);
                z2_diss1_areas(:,ar,dy) = nanmean(deoxy3D(1).sub2(1:length_diss1,cell2mat(areas(ar,row2)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(deoxy3D(2).sub1(1:length_diss2,cell2mat(areas(ar,row1)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(deoxy3D(2).sub2(1:length_diss2,cell2mat(areas(ar,row2)),dy),2);
        
                missN_s1(dy,ar,1)=numel(find(isnan(deoxy3D(1).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(deoxy3D(2).sub1(1,cell2mat(areas(ar,row1)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(deoxy3D(1).sub2(1,cell2mat(areas(ar,row2)),dy))));
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
    
    % Dyadic correlation per Conversation
    for dyad=1:numdyads
        for area=1:numareas
            a = z1_diss1_areas(:,area,dyad);
            b = z2_diss1_areas(:,area,dyad);
            [r_d1_areas(dyad,area),p_d1_areas(dyad,area)] = corr(a,b);

            c = z1_diss2_areas(:,area,dyad);
            d = z2_diss2_areas(:,area,dyad);
            [r_d2_areas(dyad,area),p_d2_areas(dyad,area)] = corr(c,d);
        end
    end     
    
    [mask_d1, ~]=fdr_bky(p_d1_areas,cutoff,'yes'); %Creates a mask for sig. areas
    [mask_d2, ~]=fdr_bky(p_d2_areas,cutoff,'yes'); %Creates a mask for sig. areas

    % Masks based on channel loss then based on FDR correction
    Sig_r_d1=r_d1_areas;
    r_FDR=r_d1_areas;
    r_FDR(~mask_d1)=NaN;
    r_FDR_ch=r_FDR;
    r_FDR_ch(~nmask(:,:,1))=NaN;
    Sig_r_d1(:,:,2)=r_FDR;
    Sig_r_d1(:,:,3)=r_FDR_ch;

    Sig_r_d2=r_d2_areas;
    r_FDR=r_d2_areas;
    r_FDR(~mask_d2)=NaN;
    r_FDR_ch=r_FDR;
    r_FDR_ch(~nmask(:,:,1))=NaN;
    Sig_r_d2(:,:,2)=r_FDR;
    Sig_r_d2(:,:,3)=r_FDR_ch;

    if writeXL
        %z-scored for later comparison
        %Sig_r_d1=atanh(Sig_r_d1);
        %Sig_r_d2=atanh(Sig_r_d2);

        load(strcat('SS_NIRS',filesep,'SS_dyads.mat'))
        dyads=table2array(dyads);
        
        if oxyOnly %file naming for type of oxy
            typeOxy='oxy_';
        else
            typeOxy='deoxy_';
        end
        
        if numareas==4 %variable names and file name for # of areas
            areas=["Subs","Dyads","mPFC","lPFC","TPJ","VIS"];
            aName='Adj_';
        elseif numareas==5
            areas=["Subs","Dyads","mPFC","lPFC","SMS","TPJ","VIS"];
            aName='OG_';
        end
        
        Sig_r_d1 = repelem(Sig_r_d1,2,1);
        Sig_r_d2 = repelem(Sig_r_d2,2,1);
        dyads=repelem(dyads,2,1);
        subs=repmat(1:2,1,numdyads)';

        Sig_r_D1=array2table([subs,dyads,Sig_r_d1(:,:,1)],'VariableNames',areas);
        r_d1_FDR=array2table([subs,dyads,Sig_r_d1(:,:,2)],'VariableNames',areas);
        r_d1_FDRch=array2table([subs,dyads,Sig_r_d1(:,:,3)],'VariableNames',areas);

        Sig_r_D2=array2table([subs,dyads,Sig_r_d2(:,:,1)],'VariableNames',areas);
        r_d2_FDR=array2table([subs,dyads,Sig_r_d2(:,:,2)],'VariableNames',areas);
        r_d2_FDRch=array2table([subs,dyads,Sig_r_d2(:,:,3)],'VariableNames',areas);
        
        if FDR
            xlName=strcat(preprocess_dir,filesep,aName,typeOxy,'D1_',num2str(cutoff),'.xls');
            writetable(Sig_r_D1,xlName,'sheet','r_vals')
            writetable(r_d1_FDR,xlName,'sheet','r_vals_FDR')
            writetable(r_d1_FDRch,xlName,'sheet','r_vals_FDRchLoss')

            xlName=strcat(preprocess_dir,filesep,aName,typeOxy,'D2_',num2str(cutoff),'.xls');
            writetable(Sig_r_D2,xlName,'sheet','r_vals')
            writetable(r_d2_FDR,xlName,'sheet','r_vals_FDR')
            writetable(r_d2_FDRch,xlName,'sheet','r_vals_FDRchLoss')
        else
            xlName=strcat(aName,typeOxy,'SS_areas.xls');
            writetable(r_d1_areas,xlName,'sheet','r_d1_areas')
            writetable(p_d1_areas,xlName,'sheet','p_d1_areas')
            writetable(r_d2_areas,xlName,'sheet','r_d2_areas')
            writetable(p_d2_areas,xlName,'sheet','p_d2_areas')
        end
    end
    save(strcat(preprocess_dir,filesep,aName,typeOxy,'SS_areas'),'Sig_r_d1','Sig_r_d2')
end
clearvars -except preprocess_dir numchans numdyads
