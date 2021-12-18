%% Instructions %%
% Dependencies for analyses: fdr_bky.m 
% Dependencies for aimaging: spm12, xjview, and mni_coordinates (in folder)
%
% Data must be preprocessed before using this script
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
oxyOnly=0;      % 0=z_totaloxy; 1=z_oxy (rarely do we look only at deoxy)
chCorr=0;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=0;     % Dyadic correlation over entire conversation per area of the brain
FDR=0;          % False discovery rate correction
writeXL=0;      % Writse the data to an excel sheet(s) in the preprocess_dir
image=0;        % To image a conversation & dyad (will prompt you for which ones)

%% Load directory & varibles to use
% You can set the directory w/ line 3 or use the get directory in line 4
%preprocess_dir= ''; %Specify path instead of selecting each run
preprocess_dir=uigetdir('','Choose Data Directory');

%% Properties
% These can be changed based on the study and what FDR cutoff is prefered
cutoff=0.01; %P-value to use as false discovery rate cut-off
length_scan=2442; %At what frame do you want to trim all subjects (shortest convo)
numdyads=54;
numchans=42;
numareas=3;

%% Compile data for easier analysis
% This will compile the two conversations of interest. Creates twelve 3D
% matrices (time,channels,dyad) based on type of oxy (deoxy,oxy,totaloxy),  
% discussion (affiliation and conflict), and subject (a and b), (3x2x2=12).
if compile
    dataprefix='0';
    %find all of the preprocessed folders
    currdir=dir(strcat(preprocess_dir,filesep,dataprefix,'*'));

% write a loop to compile the data
for i=1:length(currdir)
    dyad=currdir(i).name; %define dyad
    
    subfolder=dir(strcat(currdir(1).folder,filesep,dyad,filesep,dyad,'*',filesep,dyad,'a_affili*')); 
    if ~isempty(subfolder)
        subfile_path=strcat(subfolder.folder,filesep,subfolder.name);
        subfiles=dir(fullfile(subfile_path,'*.mat'));
        load(strcat(subfile_path,filesep,subfiles(2).name))
        length_convo=length(z_oxy);
        z_deoxy1_affil_all(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
        z_oxy1_affil_all(1:length_convo,:,i)=z_oxy(1:length_convo,:);
        z_totaloxy1_affil_all(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);
    end
    subfolder=dir(strcat(currdir(1).folder,filesep,dyad,filesep,dyad,'*',filesep,dyad,'b_affili*'));
    if ~isempty(subfolder)
        subfile_path=strcat(subfolder.folder,filesep,subfolder.name);
        subfiles=dir(fullfile(subfile_path,'*.mat'));
        load(strcat(subfile_path,filesep,subfiles(2).name))
        length_convo=length(z_oxy);
        z_deoxy2_affil_all(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
        z_oxy2_affil_all(1:length_convo,:,i)=z_oxy(1:length_convo,:);
        z_totaloxy2_affil_all(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);
    end
    
    subfolder=dir(strcat(currdir(1).folder,filesep,dyad,filesep,dyad,'*',filesep,dyad,'a_con*')); % find names of subfolder
    if ~isempty(subfolder)
        subfile_path=strcat(subfolder.folder,filesep,subfolder.name);
        subfiles=dir(fullfile(subfile_path,'*.mat'));
        load(strcat(subfile_path,filesep,subfiles(2).name))
        length_convo=length(z_oxy);
        z_deoxy1_con_all(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
        z_oxy1_con_all(1:length_convo,:,i)=z_oxy(1:length_convo,:);
        z_totaloxy1_con_all(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);
    end
    subfolder=dir(strcat(currdir(1).folder,filesep,dyad,filesep,dyad,'*',filesep,dyad,'b_con*'));
    if ~isempty(subfolder)
        subfile_path=strcat(subfolder.folder,filesep,subfolder.name);
        subfiles=dir(fullfile(subfile_path,'*.mat'));
        load(strcat(subfile_path,filesep,subfiles(2).name))
        length_convo=length(z_oxy);
        z_deoxy2_con_all(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
        z_oxy2_con_all(1:length_convo,:,i)=z_oxy(1:length_convo,:);
        z_totaloxy2_con_all(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);
    end
end

    save(strcat(preprocess_dir,filesep,'Conflict_compiled'),'z_deoxy1_affil_all','z_deoxy2_affil_all',...
'z_oxy1_affil_all','z_oxy2_affil_all','z_totaloxy1_affil_all','z_totaloxy2_affil_all','z_deoxy1_con_all','z_deoxy2_con_all',...
'z_oxy1_con_all','z_oxy2_con_all','z_totaloxy1_con_all','z_totaloxy2_con_all');

    clearvars -except preprocess_dir numdyads numchans length_scan oxyOnly chCorr areaCorr cutoff FDR writeXL image
end

%% Dyadic channel correlations 
if chCorr    
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_oxy1_con_all','z_oxy2_con_all');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = z_totaloxy1_affil_all(1:length_scan,channel,dyad);
                b = z_totaloxy2_affil_all(1:length_scan,channel,dyad);
                [r_values_affil(dyad,channel),p_values_affil(dyad,channel)] = corr(a,b);

                c = z_totaloxy1_con_all(1:length_scan,channel,dyad);
                d = z_totaloxy2_con_all(1:length_scan,channel,dyad);
                [r_values_con(dyad,channel),p_values_con(dyad,channel)] = corr(c,d);
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_totaloxy1_con_all','z_totaloxy2_con_all');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = z_oxy1_affil_all(1:length_scan,channel,dyad);
                b = z_oxy2_affil_all(1:length_scan,channel,dyad);
                [r_values_affil(dyad,channel),p_values_affil(dyad,channel)] = corr(a,b);

                c = z_totaloxy1_con_all(1:length_scan,channel,dyad);
                d = z_totaloxy2_con_all(1:length_scan,channel,dyad);
                [r_values_con(dyad,channel),p_values_con(dyad,channel)] = corr(c,d);
            end
        end
    end
    
    if FDR
        [mask_affil, ~]=fdr_bky(p_values_affil,cutoff,'yes'); 
        [mask_con, ~]=fdr_bky(p_values_con,cutoff,'yes');

        r_mask_affil=r_values_affil;
        r_mask_con=r_values_con;
        for ch=1:numchans
            r_mask_affil(~mask_affil(:,ch),ch)=NaN;
            r_mask_con(~mask_con(:,ch),ch)=NaN;
        end
        save(strcat(preprocess_dir,filesep,'CF_FDR_chCorrs'),'r_mask_affil','r_mask_con','r_values_affil','r_values_con',...
        'p_values_affil','p_values_con','mask_con','mask_affil')
    else
        save(strcat(preprocess_dir,filesep,'CF_rVals_dyadChs'),'r_values_affil','r_values_con',...
        'p_values_affil','p_values_con')
    end
    
    clearvars -except preprocess_dir numdyads numchans length_scan oxyOnly areaCorr cutoff FDR writeXL image
end

%% Mean Timecourse Synchrony per ROI per Dyad %%
if areaCorr   
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_oxy1_con_all','z_oxy2_con_all')
        
            %Get the mean for each area of interest
        z_con1_areas(:,1,:) = nanmean(z_oxy1_con_all(1:length_scan,7:14,:),2); %mpfc    
        z_con1_areas(:,2,:) = nanmean(z_oxy1_con_all(1:length_scan,[1:6,15:20],:),2); %lpfc
        z_con1_areas(:,3,:) = nanmean(z_oxy1_con_all(1:length_scan,[25:30,36:41],:),2); %tpj
        z_con2_areas(:,1,:) = nanmean(z_oxy2_con_all(1:length_scan,7:14,:),2); 
        z_con2_areas(:,2,:) = nanmean(z_oxy2_con_all(1:length_scan,[1:6,15:20],:),2);
        z_con2_areas(:,3,:) = nanmean(z_oxy2_con_all(1:length_scan,[25:30,36:41],:),2);
        %Number of missing channels per subject
        for sub=1:2
            if sub==1
                scan=z_oxy1_con_all;
            else
                scan=z_oxy2_con_all;
            end
            n_con_areas(:,1,sub) = sum(sum(isnan(scan(1:length_scan,7:14,:))))/length_scan;
            n_con_areas(:,2,sub) = sum(sum(isnan(scan(1:length_scan,[1:6,15:20],:))))/length_scan;
            n_con_areas(:,3,sub) = sum(sum(isnan(scan(1:length_scan,[25:30,36:41],:))))/length_scan;
        end
        
    else
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_totaloxy1_con_all','z_totaloxy2_con_all')
        
        %Get the mean for each area of interest
        z_con1_areas(:,1,:) = nanmean(z_totaloxy1_con_all(1:length_scan,7:14,:),2);
        z_con1_areas(:,2,:) = nanmean(z_totaloxy1_con_all(1:length_scan,[1:6,15:20],:),2);
        z_con1_areas(:,3,:) = nanmean(z_totaloxy1_con_all(1:length_scan,[25:30,36:41],:),2);
        z_con2_areas(:,1,:) = nanmean(z_totaloxy2_con_all(1:length_scan,7:14,:),2); 
        z_con2_areas(:,2,:) = nanmean(z_totaloxy2_con_all(1:length_scan,[1:6,15:20],:),2);
        z_con2_areas(:,3,:) = nanmean(z_totaloxy2_con_all(1:length_scan,[25:30,36:41],:),2);
        %Number of missing channels per subject
        for sub=1:2
            if sub==1
                scan=z_totaloxy1_con_all;
            else
                scan=z_totaloxy2_con_all;
            end
            n_con_areas(:,1,sub) = sum(sum(isnan(scan(1:length_scan,7:14,:))))/length_scan;
            n_con_areas(:,2,sub) = sum(sum(isnan(scan(1:length_scan,[1:6,15:20],:))))/length_scan;
            n_con_areas(:,3,sub) = sum(sum(isnan(scan(1:length_scan,[25:30,36:41],:))))/length_scan;
        end
    end

    %Creates mask for missing channels based on both subjects within each
    %dyad. Anything over 1/4 missing channels mask as 0. Number of channels:
    %mpfc[8/4=2chs],lpfc[12/4=3chs],tpj[12/4=3chs]
    for dyad=1:numdyads
        for area=1:numareas
            if area==1
                if n_con_areas(dyad,area,1)>1 || n_con_areas(dyad,area,2)>1
                    nmask(dyad,area)=1;
                else
                    nmask(dyad,area)=0;
                end                
            else
                if n_con_areas(dyad,area,1)>4 || n_con_areas(dyad,area,2)>4
                    nmask(dyad,area)=1;
                else
                    nmask(dyad,area)=0;
                end
            end
        end
    end
    nmask=logical(nmask);
    
    % Correlation per Conversation per ROI per Dyad 
    for area=1:numareas
        for dyad=1:numdyads
            a = z_con1_areas(:,area,dyad);
            b = z_con2_areas(:,area,dyad);
            [r_con_areas(dyad,area),p_con_areas(dyad,area)] = corr(a,b);
        end    
    end

    [mask_con, ~]=fdr_bky(p_con_areas,cutoff,'yes'); %Creates a mask for sig. areas

    Sig_r_conflict=r_con_areas;
    Sig_r_conflict(nmask)=NaN;
    Sig_r_con=Sig_r_conflict;
    Sig_r_con(~mask_con)=NaN;
    Sig_r_conflict(:,:,2)=Sig_r_con;
    Sig_r_conflict(:,:,3)=r_con_areas;

    save(strcat(preprocess_dir,filesep,'CF_areas'),'Sig_r_conflict')
    
    if writeXL
        load(strcat('CF_NIRS',filesep,'CF_dyads.mat'))
        areas=["Dyads","mPFC","lPFC","tpj"];
        Sig_r_con1=array2table([dyads,Sig_r_conflict(:,:,1)],'VariableNames',areas);
        Sig_r_con2=array2table([dyads,Sig_r_conflict(:,:,2)],'VariableNames',areas);
        r_values_con=array2table([dyads,Sig_r_conflict(:,:,3)],'VariableNames',areas);

        %File name based on FDR used
        xlName=strcat(preprocess_dir,filesep,'Sig_mask_',num2str(cutoff),'.xls'); 
        writetable(r_values_con,xlName,'sheet','r_vals')  %r-vals no mask or FDR
        writetable(Sig_r_con1,xlName,'sheet','r_vals_lostchs')  %r-vals lost chans masked
        writetable(Sig_r_con2,xlName,'sheet','r_vals_mask') %r-vals FDR & mask
    end

    clearvars -except preprocess_dir numchans numdyads image
end

%% Imaging correlated areas
if image
    load(strcat(pwd,filesep,'CF_NIRS',filesep,'CF_mnicoords.mat');
    load(strcat(preprocess_dir,filesep,'CF_FDR_chCorrs.mat'));

    %Choose the dyad and conversation you wish to image
    convoQ = 'Do you want to image the affiliation or conflict Conversation? Enter 1 for affiliation or 2 for Conflict. \n';
    convo = input(convoQ);
    
    if convo==1
        convo2img=r_mask_affil_TO;
        conversation='affil';
    else
        convo2img=r_mask_con_TO;
        conversation='con';
    end
    
    dyadQ = 'What dyad would you like to image (1 to 54)? \n';
    dyad2img = input(dyadQ);

    imgName=strcat('dyad_',num2str(dyad2img),'_',conversation,'.img');

    % Which conversation correlations you want to visualize
    convo_mask=convo2img(dyad2img,:)';

    addpath(genpath(strcat(pwd,filesep,'xjview')))
    addpath(genpath(strcat(pwd,filesep,'spm12')))
    % Make sure to change the name of the file
    nirs2img(imgName, CF_mnicoords, convo_mask, 0,0,0)
end
clear
