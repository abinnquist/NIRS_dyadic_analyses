%% Instructions %%
% Dependencies for analyses: fdr_bky.m 
% Dependencies for aimaging: spm12, xjview, and mni_coordinates (in folder)
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
areaCorr=0;     % Dyadic correlation over entire conversation per area of the brain
FDR=1;          % False discovery rate correction
writeXL=1;      % Write the data to an excel sheet(s) in the preprocess_dir
image=1;        % Will prompt you in command window for mni, conversation & dyad 

%% Load directory & varibles to use
% You can set the directory w/ line 3 or use the get directory in line 4
%preprocess_dir= ''; %Specify path instead of selecting each run
preprocess_dir=uigetdir('','Choose Data Directory');

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
    [deoxy3D,oxy3D]= compileNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans)

    %Save a .mat file to the preprocessing directory 
    save(strcat(preprocess_dir,filesep,'compiled_data'),'deoxy3D','oxy3D');

    clearvars -except preprocess_dir numdyads numchans numareas length_scan oxyOnly chCorr areaCorr cutoff FDR writeXL image
end

%% Dyadic channel correlations 
if chCorr    
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'oxy3D');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = oxy3D(1).sub1(1:length_scan,channel,dyad);
                b = oxy3D(1).sub2(1:length_scan,channel,dyad);
                [r_values_affil(dyad,channel),p_values_affil(dyad,channel)] = corr(a,b);

                c = oxy3D(2).sub1(1:length_scan,channel,dyad);
                d = oxy3D(2).sub2(1:length_scan,channel,dyad);
                [r_values_con(dyad,channel),p_values_con(dyad,channel)] = corr(c,d);
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'deoxy3D');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = deoxy3D(1).sub1(1:length_scan,channel,dyad);
                b = deoxy3D(1).sub2(1:length_scan,channel,dyad);
                [r_values_affil(dyad,channel),p_values_affil(dyad,channel)] = corr(a,b);

                c = deoxy3D(2).sub1(1:length_scan,channel,dyad);
                d = deoxy3D(2).sub2(1:length_scan,channel,dyad);
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
        
        %Get the mean for each area of interest
        for ar=1:numareas
            z_con1_areas(:,ar,:) = nanmean(oxy3D(2).sub1(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_con2_areas(:,ar,:) = nanmean(oxy3D(2).sub2(1:length_scan,cell2mat(areas(ar)),:),2); 
        end

        %Number of missing channels per subject
        for sub=1:2
            if sub==1
                scan=oxy3D(2).sub1;
            else
                scan=oxy3D(2).sub2;
            end
               %loop to count missing channels for each areas
            for nar=1:numareas
                n_con_areas(1:numdyads,nar,sub) = sum(sum(isnan(scan(1:length_scan,cell2mat(areas(nar)),:))))/length_scan;
            end
        end
            
        for dyad=1:numdyads
            for area=1:numareas
                if numareas==3
                    if area==2 || area==3
                        if n_con_areas(dyad,area,1)>2 || n_con_areas(dyad,area,2)>2
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end                
                    else
                        if n_con_areas(dyad,area,1)>1 || n_con_areas(dyad,area,2)>1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end
                    end
                elseif numareas==5
                    if area==1 || area==2
                        if n_con_areas(dyad,area,1)>=1 || n_con_areas(dyad,area,2)>=1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end                
                    else
                        if n_con_areas(dyad,area,1)>1 || n_con_areas(dyad,area,2)>1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end
                    end
                elseif numareas==6
                    if area==1 || area==5
                        if n_con_areas(dyad,area,1)>=1 || n_con_areas(dyad,area,2)>=1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end                
                    else
                        if n_con_areas(dyad,area,1)>1 || n_con_areas(dyad,area,2)>1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end
                    end
                end
            end
        end
        nmask=logical(nmask);  
    else
        load(strcat(preprocess_dir,filesep,'compiled_data.mat'),'deoxy3D')
        
        %Get the mean for each area of interest
        for ar=1:numareas
            z_con1_areas(:,ar,:) = nanmean(deoxy3D(2).sub1(1:length_scan,cell2mat(areas(ar)),:),2); 
            z_con2_areas(:,ar,:) = nanmean(deoxy3D(2).sub2(1:length_scan,cell2mat(areas(ar)),:),2); 
        end

        %Number of missing channels per subject
        for sub=1:2
            if sub==1
                scan=deoxy3D(2).sub1;
            else
                scan=deoxy3D(2).sub2;
            end
               %loop to count missing channels for each areas
            for nar=1:numareas
                n_con_areas(1:numdyads,nar,sub) = sum(sum(isnan(scan(1:length_scan,cell2mat(areas(nar)),:))))/length_scan;
            end
        end
            
        for dyad=1:numdyads
            for area=1:numareas
                if numareas==3
                    if area==2 || area==3
                        if n_con_areas(dyad,area,1)>2 || n_con_areas(dyad,area,2)>2
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end                
                    else
                        if n_con_areas(dyad,area,1)>1 || n_con_areas(dyad,area,2)>1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end
                    end
                elseif numareas==5
                    if area==1 || area==2
                        if n_con_areas(dyad,area,1)>=1 || n_con_areas(dyad,area,2)>=1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end                
                    else
                        if n_con_areas(dyad,area,1)>1 || n_con_areas(dyad,area,2)>1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end
                    end
                elseif numareas==6
                    if area==1 || area==5
                        if n_con_areas(dyad,area,1)>=1 || n_con_areas(dyad,area,2)>=1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end                
                    else
                        if n_con_areas(dyad,area,1)>1 || n_con_areas(dyad,area,2)>1
                            nmask(dyad,area)=1;
                        else
                            nmask(dyad,area)=0;
                        end
                    end
                end
            end
        end
        nmask=logical(nmask); 
    end
     
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
        
        Sig_r_con1=array2table([dyads,Sig_r_conflict(:,:,1)],'VariableNames',areas);
        Sig_r_con2=array2table([dyads,Sig_r_conflict(:,:,2)],'VariableNames',areas);
        r_values_con=array2table([dyads,Sig_r_conflict(:,:,3)],'VariableNames',areas);
        
        %File name based on FDR used
        xlName=strcat(preprocess_dir,filesep,aName,typeOxy,'con_mask_',num2str(cutoff),'.xls'); 
        writetable(r_values_con,xlName,'sheet','r_vals')  %r-vals no mask or FDR
        writetable(Sig_r_con1,xlName,'sheet','r_vals_lostchs')  %r-vals lost chans masked
        writetable(Sig_r_con2,xlName,'sheet','r_vals_mask') %r-vals FDR & mask
    end
    save(strcat(preprocess_dir,filesep,aName,typeOxy,'CF_areas'),'Sig_r_conflict')
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
