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
chCorr=1;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=1;     % Dyadic correlation over entire conversation per area of the brain
FDR=1;          % False discovery rate correction
writeXL=1;      % Write the data to an excel sheet(s) in the preprocess_dir
image=0;        % Will prompt you in command window for mni, conversation & dyad 

%% Load directory & varibles to use
% You can set the directory w/ line 3 or use the get directory in line 4
%preprocess_dir= ''; %Specify path instead of selecting each run
preprocess_dir=uigetdir('','Choose Data Directory');
addpath(genpath('helperFuncs'))

%% Properties
% These can be changed based on the study and what FDR cutoff is prefered
dataprefix='0';
ch_reject=3; % 1=no channel rejection, 2=reject noisy, 3=reject noisy & uncertain
cutoff=0.01; % Cut-off p-value to use for FDR (false discovery rate) 
length_scan=2442; % Frame to trim all subjects (shortest convo)
numdyads=54; % Number of dyads
numchans=42; % Number of channels
numareas=5;  % 3=all ch's mPFC, lPFC & TPJ, 5=subsets, 6=Lateralization

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

%% Dyadic channel correlations 
if chCorr    
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'oxy3D');
        r_values_affil=nan(numchans,numchans,numdyads);
        r_values_con=nan(numchans,numchans,numdyads);
        p_values_affil=nan(numchans,numchans,numdyads);
        p_values_con=nan(numchans,numchans,numdyads);
        
        for dyad=1:numdyads
            for chan1=1:numchans
                for chan2=1:numchans
                    a = oxy3D(1).sub1(1:length_scan,chan1,dyad);
                    b = oxy3D(1).sub2(1:length_scan,chan2,dyad);
                    [r_values_affil(chan1,chan2,dyad),p_values_affil(chan1,chan2,dyad)] = corr(a,b);
    
                    c = oxy3D(2).sub1(1:length_scan,chan1,dyad);
                    d = oxy3D(2).sub2(1:length_scan,chan2,dyad);
                    [r_values_con(chan1,chan2,dyad),p_values_con(chan1,chan2,dyad)] = corr(c,d);
                end
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'deoxy3D');
        r_values_affil=nan(numchans,numchans,numdyads);
        r_values_con=nan(numchans,numchans,numdyads);
        p_values_affil=nan(numchans,numchans,numdyads);
        p_values_con=nan(numchans,numchans,numdyads);
        
        for dyad=1:numdyads
            for chan1=1:numchans
                for chan2=1:numchans
                    a = deoxy3D(1).sub1(1:length_scan,chan1,dyad);
                    b = deoxy3D(1).sub2(1:length_scan,chan2,dyad);
                    [r_values_affil(chan1,chan2,dyad),p_values_affil(chan1,chan2,dyad)] = corr(a,b);
    
                    c = deoxy3D(2).sub1(1:length_scan,chan1,dyad);
                    d = deoxy3D(2).sub2(1:length_scan,chan2,dyad);
                    [r_values_con(chan1,chan2,dyad),p_values_con(chan1,chan2,dyad)] = corr(c,d);
                end
            end
        end
    end

    %for paired channels corrs only
    for dy=1:numdyads
        for ch=1:numchans
            r_paired_affil(ch,dy)=r_values_affil(ch,ch,dy);
            r_paired_con(ch,dy)=r_values_con(ch,ch,dy);
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
        'p_values_affil','p_values_con','mask_con','mask_affil','r_paired_affil','r_paired_con')
    else
        save(strcat(preprocess_dir,filesep,'CF_rVals_dyadChs'),'r_values_affil','r_values_con',...
        'p_values_affil','p_values_con','r_paired_affil','r_paired_con')
    end
    
    clearvars -except preprocess_dir numdyads numchans numareas length_scan oxyOnly areaCorr cutoff FDR writeXL image
end

%% Mean Timecourse Synchrony per ROI per Dyad %%
if areaCorr
    if numareas==3 
        areas={7:14;[1:6,15:20];[25:30,36:41]}; %mPFC, lPFC, TPJ
    elseif numareas==5
        areas={[11,13];[8,10];[4,6,15,19];[1,3,18,20];[25,28,37,39]}; %vmPFC, dmPFC, vlPFC, dlPFC, TPJ
    elseif numareas==6
        areas={[8,11];1:6;25:29;[10,13];15:20;36:40};%L-mPFC, L-lPFC, L-TPJ, R-mPFC, R-lPFC, R-TPJ
    end  
   
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'oxy3D')

        %Get the mean for each area of interest, for every timepoint
        for ar=1:numareas
            z_aff1_areas(:,ar,:) = nanmean(oxy3D(1).sub1(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_aff2_areas(:,ar,:) = nanmean(oxy3D(1).sub2(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_con1_areas(:,ar,:) = nanmean(oxy3D(2).sub1(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_con2_areas(:,ar,:) = nanmean(oxy3D(2).sub2(1:length_scan,cell2mat(areas(ar)),:),2); 
        end

        %Get missing channels per area of interest
        for convo=1:2
            for dy=1:numdyads
                for ar=1:numareas
                    missN_s1(dy,ar,convo) = numel(find(isnan(oxy3D(convo).sub1(1,cell2mat(areas(ar)),dy)))); 
                    missN_s2(dy,ar,convo) = numel(find(isnan(oxy3D(convo).sub2(1,cell2mat(areas(ar)),dy)))); 
                end
            end
        end
        nmask = cfMissingch(missN_s1,missN_s2,numdyads,numareas);
        nmask=logical(nmask);  

    else
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'deoxy3D')
        
        %Get the mean for each area of interest
        for ar=1:numareas
            z_aff1_areas(:,ar,:) = nanmean(deoxy3D(1).sub1(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_aff2_areas(:,ar,:) = nanmean(deoxy3D(1).sub2(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_con1_areas(:,ar,:) = nanmean(deoxy3D(2).sub1(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_con2_areas(:,ar,:) = nanmean(deoxy3D(2).sub2(1:length_scan,cell2mat(areas(ar)),:),2); 
        end

        %Number of missing channels per subject
        for convo=1:2
            for dy=1:numdyads
                for ar=1:numareas
                    missN_s1(dy,ar,convo) = numel(find(isnan(deoxy3D(convo).sub1(1,cell2mat(areas(ar)),dy)))); 
                    missN_s2(dy,ar,convo) = numel(find(isnan(deoxy3D(convo).sub2(1,cell2mat(areas(ar)),dy)))); 
                end
            end
        end
        nmask = cfMissingch(missN_s1,missN_s2,numdyads,numareas);
        nmask=logical(nmask);   
    end
     
    % Correlation of ROI x dyad x convo
    for area=1:numareas
        for dyad=1:numdyads
            a = z_aff1_areas(:,area,dyad);
            b = z_aff2_areas(:,area,dyad);
            [r_aff_areas(dyad,area),p_aff_areas(dyad,area)] = corr(a,b);

            c = z_con1_areas(:,area,dyad);
            d = z_con2_areas(:,area,dyad);
            [r_con_areas(dyad,area),p_con_areas(dyad,area)] = corr(c,d);
        end    
    end

    [mask_aff, ~]=fdr_bky(p_aff_areas,cutoff,'yes'); %Creates a mask for sig. areas
    [mask_con, ~]=fdr_bky(p_con_areas,cutoff,'yes'); %Creates a mask for sig. areas

    Sig_r_affiliation=r_aff_areas;
    Sig_r_affiliation(nmask(:,:,1))=NaN;
    Sig_r_aff=Sig_r_affiliation;
    Sig_r_aff(~mask_aff)=NaN;
    Sig_r_affiliation(:,:,2)=Sig_r_aff;
    Sig_r_affiliation(:,:,3)=r_aff_areas;

    Sig_r_conflict=r_con_areas;
    Sig_r_conflict(nmask(:,:,2))=NaN;
    Sig_r_con=Sig_r_conflict;
    Sig_r_con(~mask_con)=NaN;
    Sig_r_conflict(:,:,2)=Sig_r_con;
    Sig_r_conflict(:,:,3)=r_con_areas;
     
    if writeXL
        Sig_r_conflict=atanh(Sig_r_conflict);%z-scored r-values
        load(strcat('CF_NIRS',filesep,'CF_dyads.mat'))
        
        if oxyOnly
            typeOxy='oxy_'; %file naming for type of oxy
        else
            typeOxy='deoxy_';
        end
            
        if numareas==3 %variable names and file name for # of areas
            areas=["Dyads","mPFC","lPFC","tpj"];
            aName='OG_';
        elseif numareas==5
            areas=["Dyads","vmPFC","dmPFC","vlPFC","dlPFC","tpj"];
            aName='VD_';
        elseif numareas==6
            areas=["Dyads","L-mPFC","L-lPFC","L-tpj","R-mPFC","R-lPFC","R-tpj"];
            aName='Lat_';
        end
        
        Sig_r_aff1=array2table([dyads,Sig_r_affiliation(:,:,1)],'VariableNames',areas);
        Sig_r_aff2=array2table([dyads,Sig_r_affiliation(:,:,2)],'VariableNames',areas);
        r_values_aff=array2table([dyads,Sig_r_affiliation(:,:,3)],'VariableNames',areas);

        Sig_r_con1=array2table([dyads,Sig_r_conflict(:,:,1)],'VariableNames',areas);
        Sig_r_con2=array2table([dyads,Sig_r_conflict(:,:,2)],'VariableNames',areas);
        r_values_con=array2table([dyads,Sig_r_conflict(:,:,3)],'VariableNames',areas);
        
        %File name based on FDR used
        xlName=strcat(preprocess_dir,filesep,aName,typeOxy,'con_mask_',num2str(cutoff),'.xls'); 
        writetable(r_values_con,xlName,'sheet','r_vals')  %r-vals no mask or FDR
        writetable(Sig_r_con1,xlName,'sheet','r_vals_lostchs')  %r-vals lost chans masked
        writetable(Sig_r_con2,xlName,'sheet','r_vals_mask') %r-vals FDR & mask
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
