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
chCorr=0;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=1;     % Dyadic correlation over entire conversation per area of the brain
FDR=1;          % False discovery rate correction
writeXL=1;      % Write the data to an excel sheet(s) in the preprocess_dir
image=0;        % Will prompt you in command window for mni, conversation & dyad 

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
numareas=6;  % 3=all ch's mPFC, lPFC & TPJ, 5=subsets, 6=Lateralization

%% Compile data for easier analysis
% This will compile the two conversations of interest. Creates twelve 3D
% matrices (time,channels,dyad) based on type of oxy (deoxy,oxy,totaloxy),  
% discussion (affiliation and conflict), and subject (a and b), (3x2x2=12).
if compile
    [z_deoxy1_1,z_oxy1_1,z_totaloxy1_1,z_deoxy1_2,z_oxy1_2,...
    z_totaloxy1_2,z_deoxy2_1,z_oxy2_1,z_totaloxy2_1,z_deoxy2_2,...
    z_oxy2_2,z_totaloxy2_2]= compileNIRSdata(preprocess_dir,dataprefix,ch_reject);

    z_deoxy1_affil=z_deoxy1_1;
    z_deoxy2_affil=z_deoxy1_2;
    z_deoxy1_con=z_deoxy2_1;
    z_deoxy2_con=z_deoxy2_2;
    z_oxy1_affil=z_oxy1_1;
    z_oxy2_affil=z_oxy1_2;
    z_oxy1_con=z_oxy2_1;
    z_oxy2_con=z_oxy2_2;
    z_totaloxy1_affil=z_totaloxy1_1;
    z_totaloxy2_affil=z_totaloxy1_2;
    z_totaloxy1_con=z_totaloxy2_1;
    z_totaloxy2_con=z_totaloxy2_2;

    save(strcat(preprocess_dir,filesep,'Conflict_compiled'),'z_deoxy1_affil','z_deoxy2_affil',...
    'z_oxy1_affil','z_oxy2_affil','z_totaloxy1_affil','z_totaloxy2_affil','z_deoxy1_con','z_deoxy2_con',...
    'z_oxy1_con','z_oxy2_con','z_totaloxy1_con','z_totaloxy2_con');

    clearvars -except preprocess_dir numdyads numchans numareas length_scan oxyOnly chCorr areaCorr cutoff FDR writeXL image
end

%% Dyadic channel correlations 
if chCorr    
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_oxy1_con','z_oxy2_con');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = z_totaloxy1_affil(1:length_scan,channel,dyad);
                b = z_totaloxy2_affil(1:length_scan,channel,dyad);
                [r_values_affil(dyad,channel),p_values_affil(dyad,channel)] = corr(a,b);

                c = z_totaloxy1_con(1:length_scan,channel,dyad);
                d = z_totaloxy2_con(1:length_scan,channel,dyad);
                [r_values_con(dyad,channel),p_values_con(dyad,channel)] = corr(c,d);
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_totaloxy1_con','z_totaloxy2_con');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = z_oxy1_affil(1:length_scan,channel,dyad);
                b = z_oxy2_affil(1:length_scan,channel,dyad);
                [r_values_affil(dyad,channel),p_values_affil(dyad,channel)] = corr(a,b);

                c = z_totaloxy1_con(1:length_scan,channel,dyad);
                d = z_totaloxy2_con(1:length_scan,channel,dyad);
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
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_oxy1_con','z_oxy2_con')
        
        if numareas==6 %lateralization of PFC & TPJ
            %Get the mean for each area of interest
            %left side
            z_con1_areas(:,1,:) = nanmean(z_oxy1_con(1:length_scan,[8,11],:),2); %L-mpfc    
            z_con1_areas(:,2,:) = nanmean(z_oxy1_con(1:length_scan,1:6,:),2); %L-lpfc
            z_con1_areas(:,3,:) = nanmean(z_oxy1_con(1:length_scan,36:40,:),2); %L-tpj
                %right side
            z_con1_areas(:,4,:) = nanmean(z_oxy1_con(1:length_scan,[10,13],:),2); %R-mpfc    
            z_con1_areas(:,5,:) = nanmean(z_oxy1_con(1:length_scan,15:20,:),2); %R-lpfc
            z_con1_areas(:,6,:) = nanmean(z_oxy1_con(1:length_scan,25:29,:),2); %R-tpj

                        %left side
            z_con2_areas(:,1,:) = nanmean(z_oxy2_con(1:length_scan,[8,11],:),2); %mpfc    
            z_con2_areas(:,2,:) = nanmean(z_oxy2_con(1:length_scan,1:6,:),2); %lpfc
            z_con2_areas(:,3,:) = nanmean(z_oxy2_con(1:length_scan,36:40,:),2); %tpj
                %right side
            z_con2_areas(:,4,:) = nanmean(z_oxy2_con(1:length_scan,[10,13],:),2); %mpfc    
            z_con2_areas(:,5,:) = nanmean(z_oxy2_con(1:length_scan,15:20,:),2); %lpfc
            z_con2_areas(:,6,:) = nanmean(z_oxy2_con(1:length_scan,25:29,:),2); %tpj
            
            %Number of missing channels per subject
            for sub=1:2
                if sub==1
                    scan=z_oxy1_con;
                else
                    scan=z_oxy2_con;
                end
                n_con_areas(:,1,sub) = sum(sum(isnan(scan(1:length_scan,[10,13],:))))/length_scan;
                n_con_areas(:,2,sub) = sum(sum(isnan(scan(1:length_scan,15:20,:))))/length_scan;
                n_con_areas(:,3,sub) = sum(sum(isnan(scan(1:length_scan,25:29,:))))/length_scan;
                n_con_areas(:,4,sub) = sum(sum(isnan(scan(1:length_scan,[8,11],:))))/length_scan;
                n_con_areas(:,5,sub) = sum(sum(isnan(scan(1:length_scan,1:6,:))))/length_scan;
                n_con_areas(:,6,sub) = sum(sum(isnan(scan(1:length_scan,36:40,:))))/length_scan;
            end
            
            for dyad=1:numdyads
                for area=1:numareas
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
                end
            end
            nmask=logical(nmask);
            
        elseif numareas==5 % Smaller subsets of PFC & TPJ
            z_con1_areas(:,1,:) = nanmean(z_oxy1_con(1:length_scan,[11,13],:),2); %vmPFC    
            z_con1_areas(:,2,:) = nanmean(z_oxy1_con(1:length_scan,[8,10],:),2); %dmPFC
            z_con1_areas(:,3,:) = nanmean(z_oxy1_con(1:length_scan,[4,6,15,19],:),2); %vlPFC
            z_con1_areas(:,4,:) = nanmean(z_oxy1_con(1:length_scan,[1,3,18,20],:),2); %dlPFC    
            z_con1_areas(:,5,:) = nanmean(z_oxy1_con(1:length_scan,[25,28,37,39],:),2); %TPJ

            z_con2_areas(:,1,:) = nanmean(z_oxy2_con(1:length_scan,[11,13],:),2); %vmPFC    
            z_con2_areas(:,2,:) = nanmean(z_oxy2_con(1:length_scan,[8,10],:),2); %dmPFC
            z_con2_areas(:,3,:) = nanmean(z_oxy2_con(1:length_scan,[4,6,15,19],:),2); %vlPFC
            z_con2_areas(:,4,:) = nanmean(z_oxy2_con(1:length_scan,[1,3,18,20],:),2); %dlPFC    
            z_con2_areas(:,5,:) = nanmean(z_oxy2_con(1:length_scan,[25,28,37,39],:),2); %TPJ
            
            for sub=1:2
                if sub==1
                    scan=z_oxy1_con;
                else
                    scan=z_oxy2_con;
                end
                n_con_areas(:,1,sub) = sum(sum(isnan(scan(1:length_scan,[11,13],:))))/length_scan;
                n_con_areas(:,2,sub) = sum(sum(isnan(scan(1:length_scan,[8,10],:))))/length_scan;
                n_con_areas(:,3,sub) = sum(sum(isnan(scan(1:length_scan,[4,6,15,19],:))))/length_scan;
                n_con_areas(:,4,sub) = sum(sum(isnan(scan(1:length_scan,[1,3,18,20],:))))/length_scan;
                n_con_areas(:,5,sub) = sum(sum(isnan(scan(1:length_scan,[25,28,37,39],:))))/length_scan;
            end
            
            for dyad=1:numdyads
                for area=1:numareas
                    if area==1
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
            nmask=logical(nmask);
          
        elseif numareas==3 %including all channels fro pfc & tpj
            z_con1_areas(:,1,:) = nanmean(z_oxy1_con(1:length_scan,7:14,:),2);
            z_con1_areas(:,2,:) = nanmean(z_oxy1_con(1:length_scan,[1:6,15:20],:),2);
            z_con1_areas(:,3,:) = nanmean(z_oxy1_con(1:length_scan,[25:30,36:41],:),2);

            z_con2_areas(:,1,:) = nanmean(z_oxy2_con(1:length_scan,7:14,:),2); 
            z_con2_areas(:,2,:) = nanmean(z_oxy2_con(1:length_scan,[1:6,15:20],:),2);
            z_con2_areas(:,3,:) = nanmean(z_oxy2_con(1:length_scan,[25:30,36:41],:),2);
            
            for sub=1:2
                if sub==1
                    scan=z_deoxy1_con;
                else
                    scan=z_deoxy2_con;
                end
                n_con_areas(:,1,sub) = sum(sum(isnan(scan(1:length_scan,7:14,:))))/length_scan;
                n_con_areas(:,2,sub) = sum(sum(isnan(scan(1:length_scan,[1:6,15:20],:))))/length_scan;
                n_con_areas(:,3,sub) = sum(sum(isnan(scan(1:length_scan,[25:30,36:41],:))))/length_scan;
            end
        end 
    else
        load(strcat(preprocess_dir,filesep,'Conflict_compiled.mat'),'z_deoxy1_con','z_deoxy2_con')
        
        %Get the mean for each area of interest
        z_con1_areas(:,1,:) = nanmean(z_deoxy1_con(1:length_scan,7:14,:),2);
        z_con1_areas(:,2,:) = nanmean(z_deoxy1_con(1:length_scan,[1:6,15:20],:),2);
        z_con1_areas(:,3,:) = nanmean(z_deoxy1_con(1:length_scan,[25:30,36:41],:),2);
        
        z_con2_areas(:,1,:) = nanmean(z_deoxy2_con(1:length_scan,7:14,:),2); 
        z_con2_areas(:,2,:) = nanmean(z_deoxy2_con(1:length_scan,[1:6,15:20],:),2);
        z_con2_areas(:,3,:) = nanmean(z_deoxy2_con(1:length_scan,[25:30,36:41],:),2);
        %Number of missing channels per subject
        for sub=1:2
            if sub==1
                scan=z_deoxy1_con;
            else
                scan=z_deoxy2_con;
            end
            n_con_areas(:,1,sub) = sum(sum(isnan(scan(1:length_scan,7:14,:))))/length_scan;
            n_con_areas(:,2,sub) = sum(sum(isnan(scan(1:length_scan,[1:6,15:20],:))))/length_scan;
            n_con_areas(:,3,sub) = sum(sum(isnan(scan(1:length_scan,[25:30,36:41],:))))/length_scan;
        end
        
       %Creates mask for missing channels 
        for dyad=1:numdyads
            for area=1:numareas
                if area==1
                    if n_con_areas(dyad,area,1)>2 || n_con_areas(dyad,area,2)>2
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

    save(strcat(preprocess_dir,filesep,'CF_areas'),'Sig_r_conflict')
     
    if writeXL
        %z-scored for later comparison
        Sig_r_conflict=atanh(Sig_r_conflict);
        
        load(strcat('CF_NIRS',filesep,'CF_dyads.mat'))
        if area==3
            areas=["Dyads","mPFC","lPFC","tpj"];
        elseif area==5
            areas=["Dyads","vmPFC","dmPFC","vlPFC","dlPFC","tpj"];
        elseif area==6
            areas=["Dyads","L-mPFC","L-lPFC","L-tpj","R-mPFC","R-lPFC","R-tpj"];
        end
        
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
