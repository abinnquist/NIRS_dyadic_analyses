%% Instructions %%
% Dependencies for analyses: fdr_bky.m 
% Dependencies for imaging: spm12, xjview, and mni_coordinates (in folder)
% 
% Data must be preprocessed before using this script (see
% https://github.com/abinnquist/fNIRSpreProcessing)
%
% Before running the script make sure to choose the analyses you want to run 
% by turning on (1), or turning off (0), the tasks you want to perform. 
% You can change the properties for a different trim time, more 
% strict/leneint FDR cutoff and/or channel rejection based on the three 
% outputs from preprocessing,
% different trim time.
%
% If you prefer to specify your data path instead of selecting each time
% comment out the uigetdir command and uncomment the command above it. Make
% sure to specify the location of the preprocessed NIRS data
%% Analyses to run
% Set to zero if you do not want to perform the task
compile=1;      % If the data has yet to be compiled into all dyad matices
oxyOnly=1;      % 0=z_deoxy; 1=z_oxy
chCorr=0;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=0;     % Dyadic correlation over entire conversation per area of the brain
FDR=0;          % False discovery rate correction
writeXL=0;      % If you want to write the data to an excel sheet(s)
image=0;        % Will prompt you in command window for mni, conversation & dyad 

%% Set the Directory
% You can set the directory w/ line 3 or use the get directory in line 4
% preprocess_dir= ''; %Make sure to specify correct location
preprocess_dir=uigetdir('','Choose Data Directory');

%% Properties
% These can be changed based on the study and what FDR cutoff is prefered
dataprefix='SS';
ch_reject=3;  % 1=no channel rejection, 2=reject noisy, 3=reject noisy & uncertain
cutoff=0.1; % Cut-off p-value to use for FDR (false discovery rate) 
length_diss1=1853; % Can be changed to the shortest conversation
length_diss2=1810; % Can be changed to the shortest conversation
numdyads=52; % Number of dyads
numchans=42; % Number of channels
numareas=6; % 4=Right Lateralization; 5=mPFC, lPFC, PMC, SMS & TPJ; 6=vmPFC, dmPFC, vlPFC, SMS, TPJ

%% Compile subject data %%
% Compiles the two conversations of interest. Creates a 3D
% matrix for type of oxy (deoxy & oxy), discussion (1 and 2),  
% and subject (1 and 2) for eight (2x2x2=12) 3D matrices (time,channels,dyad)
if compile
    numScans=2;
    [deoxy3D,oxy3D]= compileNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans);

    %Save a .mat file to the preprocessing directory 
    save(strcat(preprocess_dir,filesep,'SS_compiled'),'deoxy3D','oxy3D');
    
    clearvars -except preprocess_dir numdyads numchans numareas length_diss1 length_diss2 oxyOnly chCorr areaCorr cutoff FDR image writeXL
end

%% Compute basic channel correlations %%
% Computes the matched channel correlations for each dyad and conversation. 
if chCorr
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'oxy3D');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = oxy3D(1).sub1(1:length_diss1,channel,dyad);
                b = oxy3D(1).sub2(1:length_diss1,channel,dyad);
                [r_values_diss1(dyad,channel),p_values_diss1(dyad,channel)] = corr(a,b);

                c = oxy3D(2).sub1(1:length_diss2,channel,dyad);
                d = oxy3D(2).sub2(1:length_diss2,channel,dyad);
                [r_values_diss2(dyad,channel),p_values_diss2(dyad,channel)] = corr(c,d);
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'deoxy3D');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = deoxy3D(1).sub1(:,channel,dyad);
                b = deoxy3D(1).sub2(:,channel,dyad);
                [r_values_diss1(dyad,channel),p_values_diss1(dyad,channel)] = corr(a,b);

                c = deoxy3D(2).sub1(:,channel,dyad);
                d = deoxy3D(2).sub2(:,channel,dyad);
                [r_values_diss2(dyad,channel),p_values_diss2(dyad,channel)] = corr(c,d);
            end
        end
    end
    %Uses the FDR_bky correction based on the cutoff specified in properties
    if FDR
        [mask_diss1, ~]=fdr_bky(p_values_diss1,cutoff,'yes'); 
        [mask_diss2, ~]=fdr_bky(p_values_diss2,cutoff,'yes');

        r_mask_diss1=r_values_diss1;
        r_mask_diss2=r_values_diss2;
        for ch=1:numchans
            r_mask_diss1(~mask_diss1(:,ch),ch)=NaN;
            r_mask_diss2(~mask_diss2(:,ch),ch)=NaN;
        end
        FDR_r_mask=r_mask_diss1;
        FDR_r_mask(:,:,2)=r_mask_diss2;
        save(strcat(preprocess_dir,filesep,'FDR_r_mask'),'FDR_r_mask')
        save(strcat(preprocess_dir,filesep,'SS_FDR_chCorrs'),'r_values_diss1','r_values_diss2',...
        'p_values_diss1','p_values_diss2','mask_diss1','mask_diss2')
    else
        save(strcat(preprocess_dir,filesep,'SS_rVals_dyadChs'),'r_values_diss1','r_values_diss2',...
        'p_values_diss1','p_values_diss2')
    end
    %If you want to output an excel file with r-values
    if writeXL
        if FDR
            writetable(r_mask_diss1,'FDR_SS_channels_0.1.xls','sheet','FDR_r_vals_diss1')
            writetable(r_mask_diss2,'FDR_SS_channels_0.1.xls','sheet','FDR_r_vals_diss2')
        else
            writetable(r_values_diss1,'Sig_SS_channels.xls','sheet','r_values_diss1')
            writetable(p_values_diss1,'Sig_SS_channels.xls','sheet','p_values_diss1')
            writetable(r_values_diss2,'Sig_SS_channels.xls','sheet','r_values_diss2')
            writetable(p_values_diss2,'Sig_SS_channels.xls','sheet','p_values_diss2')
        end
    end

    clearvars -except preprocess_dir numdyads numchans numareas length_diss1 length_diss2 oxyOnly chCorr areaCorr cutoff FDR image writeXL
end

%% Mean Timecourse Synchrony per ROI per Dyad %%
% Averages activiation in area specific channels then computes dyadic 
% correlation of those areas of the brain 
if areaCorr           
    if numareas==4 %Right lateralization
        area1=[2,8,12]; %medial prefrontal cortex
        area2=[7,30:31];  %lateral prefrontal cortex
        area3=[29,33:35]; %somatosensory cortex
        area4=[36:40,42]; %temporoparietal junction
    elseif numareas==5 %areas including all channels
        area1=[1:3,5,8,10:12]; %medial prefrontal cortex
        area2=[4,6,7,9,30:31];  %lateral prefrontal cortex
        area3=[13,24:25,32]; %premotor cortex
        area4=[14:16,26:29,33:35]; %%somatosensory cortex
        area5=[17:21,23,36:40,42]; %temporoparietal junction
    elseif numareas==6 %Sub-areas with reliable channels
        area1=1:3; %ventral medialPFC
        area2=11:12;  %dorsal medial PFC
        area3=[4,7,30:31]; %ventral lateral PFC
        area4=[14:16,33:35]; %somatosensory cortex
        area5=[17:19,36:38]; %anterior temporoparietal junction
        area6=[21,23,40,42]; %posterior TPJ
    end   
    
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'oxy3D');

        z1_diss1_areas(:,1,:) = nanmean(oxy3D(1).sub1(1:length_diss1,area1,:),2);   
        z1_diss1_areas(:,2,:) = nanmean(oxy3D(1).sub1(1:length_diss1,area2,:),2); 
        z1_diss1_areas(:,3,:) = nanmean(oxy3D(1).sub1(1:length_diss1,area3,:),2); 
        z1_diss1_areas(:,4,:) = nanmean(oxy3D(1).sub1(1:length_diss1,area4,:),2); 

        z2_diss1_areas(:,1,:) = nanmean(oxy3D(1).sub2(1:length_diss1,area1,:),2); 
        z2_diss1_areas(:,2,:) = nanmean(oxy3D(1).sub2(1:length_diss1,area2,:),2);
        z2_diss1_areas(:,3,:) = nanmean(oxy3D(1).sub2(1:length_diss1,area3,:),2);
        z2_diss1_areas(:,4,:) = nanmean(oxy3D(1).sub2(1:length_diss1,area4,:),2);

        z1_diss2_areas(:,1,:) = nanmean(oxy3D(2).sub1(1:length_diss2,area1,:),2);  
        z1_diss2_areas(:,2,:) = nanmean(oxy3D(2).sub1(1:length_diss2,area2,:),2);
        z1_diss2_areas(:,3,:) = nanmean(oxy3D(2).sub1(1:length_diss2,area3,:),2); 
        z1_diss2_areas(:,4,:) = nanmean(oxy3D(2).sub1(1:length_diss2,area4,:),2);

        z2_diss2_areas(:,1,:) = nanmean(oxy3D(2).sub2(1:length_diss2,area1,:),2); 
        z2_diss2_areas(:,2,:) = nanmean(oxy3D(2).sub2(1:length_diss2,area2,:),2);
        z2_diss2_areas(:,3,:) = nanmean(oxy3D(2).sub2(1:length_diss2,area3,:),2);
        z2_diss2_areas(:,4,:) = nanmean(oxy3D(2).sub2(1:length_diss2,area4,:),2);
        
        if numareas==5
            z1_diss1_areas(:,5,:) = nanmean(oxy3D(1).sub1(1:length_diss1,area5,:),2);
            z2_diss1_areas(:,5,:) = nanmean(oxy3D(1).sub2(1:length_diss1,area5,:),2);
            z1_diss2_areas(:,5,:) = nanmean(oxy3D(2).sub1(1:length_diss2,area5,:),2);
            z2_diss2_areas(:,5,:) = nanmean(oxy3D(2).sub2(1:length_diss2,area5,:),2);
        elseif numareas==6
            z1_diss1_areas(:,5,:) = nanmean(oxy3D(1).sub1(1:length_diss1,area5,:),2);
            z2_diss1_areas(:,5,:) = nanmean(oxy3D(1).sub2(1:length_diss1,area5,:),2);
            z1_diss2_areas(:,5,:) = nanmean(oxy3D(2).sub1(1:length_diss2,area5,:),2);
            z2_diss2_areas(:,5,:) = nanmean(oxy3D(2).sub2(1:length_diss2,area5,:),2);
            
            z1_diss1_areas(:,6,:) = nanmean(oxy3D(1).sub1(1:length_diss1,area6,:),2);
            z2_diss1_areas(:,6,:) = nanmean(oxy3D(1).sub2(1:length_diss1,area6,:),2);
            z1_diss2_areas(:,6,:) = nanmean(oxy3D(2).sub1(1:length_diss2,area6,:),2);
            z2_diss2_areas(:,6,:) = nanmean(oxy3D(2).sub2(1:length_diss2,area6,:),2);
        end

        %Number of missing channels per subject, per area of interest
        for dc=1:2
            if dc==1
                for sub=1:2
                    if sub==1
                        scan=oxy3D(1).sub1;
                    else
                        scan=oxy3D(1).sub2;
                    end
                    n_m1_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss1,area1,:))))/length_diss1;
                    n_m1_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss1,area2,:))))/length_diss1;
                    n_m1_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss1,area3,:))))/length_diss1;
                    n_m1_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss1,area4,:))))/length_diss1; 
                    if numareas==5
                        n_m1_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,area5,:))))/length_diss1; 
                    elseif numareas==6
                        n_m1_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,area5,:))))/length_diss1; 
                        n_m1_areas(:,6,sub) = sum(sum(isnan(scan(1:length_diss1,area6,:))))/length_diss1; 
                    end
                end
            else
                for sub=1:2
                    if sub==1
                        scan=oxy3D(2).sub1;
                    else
                        scan=oxy3D(2).sub2;
                    end
                    n_m2_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss2,area1,:))))/length_diss2;
                    n_m2_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss2,area2,:))))/length_diss2;
                    n_m2_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss2,area3,:))))/length_diss2;
                    n_m2_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss2,area4,:))))/length_diss2; 
                    if numareas==5
                        n_m2_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss2,area5,:))))/length_diss2; 
                    elseif numareas==6
                        n_m2_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss2,area5,:))))/length_diss2; 
                        n_m2_areas(:,6,sub) = sum(sum(isnan(scan(1:length_diss2,area6,:))))/length_diss2; 
                    end
                end
            end
        end
      
    else
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'deoxy3D');

        z1_diss1_areas(:,1,:) = nanmean(deoxy3D(1).sub1(1:length_diss1,area1,:),2);   
        z1_diss1_areas(:,2,:) = nanmean(deoxy3D(1).sub1(1:length_diss1,area2,:),2); 
        z1_diss1_areas(:,3,:) = nanmean(deoxy3D(1).sub1(1:length_diss1,area3,:),2); 
        z1_diss1_areas(:,4,:) = nanmean(deoxy3D(1).sub1(1:length_diss1,area4,:),2); 

        z2_diss1_areas(:,1,:) = nanmean(deoxy3D(1).sub2(1:length_diss1,area1,:),2); 
        z2_diss1_areas(:,2,:) = nanmean(deoxy3D(1).sub2(1:length_diss1,area2,:),2);
        z2_diss1_areas(:,3,:) = nanmean(deoxy3D(1).sub2(1:length_diss1,area3,:),2);
        z2_diss1_areas(:,4,:) = nanmean(deoxy3D(1).sub2(1:length_diss1,area4,:),2);

        z1_diss2_areas(:,1,:) = nanmean(deoxy3D(2).sub1(1:length_diss2,area1,:),2);  
        z1_diss2_areas(:,2,:) = nanmean(deoxy3D(2).sub1(1:length_diss2,area2,:),2);
        z1_diss2_areas(:,3,:) = nanmean(deoxy3D(2).sub1(1:length_diss2,area3,:),2); 
        z1_diss2_areas(:,4,:) = nanmean(deoxy3D(2).sub1(1:length_diss2,area4,:),2);

        z2_diss2_areas(:,1,:) = nanmean(deoxy3D(2).sub2(1:length_diss2,area1,:),2); 
        z2_diss2_areas(:,2,:) = nanmean(deoxy3D(2).sub2(1:length_diss2,area2,:),2);
        z2_diss2_areas(:,3,:) = nanmean(deoxy3D(2).sub2(1:length_diss2,area3,:),2);
        z2_diss2_areas(:,4,:) = nanmean(deoxy3D(2).sub2(1:length_diss2,area4,:),2);
        
        if numareas==5
            z1_diss1_areas(:,5,:) = nanmean(deoxy3D(1).sub1(1:length_diss1,area5,:),2);
            z2_diss1_areas(:,5,:) = nanmean(deoxy3D(1).sub2(1:length_diss1,area5,:),2);
            z1_diss2_areas(:,5,:) = nanmean(deoxy3D(2).sub1(1:length_diss2,area5,:),2);
            z2_diss2_areas(:,5,:) = nanmean(deoxy3D(2).sub2(1:length_diss2,area5,:),2);
        elseif numareas==6
            z1_diss1_areas(:,5,:) = nanmean(deoxy3D(1).sub1(1:length_diss1,area5,:),2);
            z2_diss1_areas(:,5,:) = nanmean(deoxy3D(1).sub2(1:length_diss1,area5,:),2);
            z1_diss2_areas(:,5,:) = nanmean(deoxy3D(2).sub1(1:length_diss2,area5,:),2);
            z2_diss2_areas(:,5,:) = nanmean(deoxy3D(2).sub2(1:length_diss2,area5,:),2);
            
            z1_diss1_areas(:,6,:) = nanmean(deoxy3D(1).sub1(1:length_diss1,area6,:),2);
            z2_diss1_areas(:,6,:) = nanmean(deoxy3D(1).sub2(1:length_diss1,area6,:),2);
            z1_diss2_areas(:,6,:) = nanmean(deoxy3D(2).sub1(1:length_diss2,area6,:),2);
            z2_diss2_areas(:,6,:) = nanmean(deoxy3D(2).sub2(1:length_diss2,area6,:),2);
        end

        %Number of missing channels per subject, per area of interest
        for dc=1:2
            if dc==1
                for sub=1:2
                    if sub==1
                        scan=deoxy3D(1).sub1;
                    else
                        scan=deoxy3D(1).sub2;
                    end
                    n_m1_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss1,area1,:))))/length_diss1;
                    n_m1_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss1,area2,:))))/length_diss1;
                    n_m1_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss1,area3,:))))/length_diss1;
                    n_m1_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss1,area4,:))))/length_diss1; 
                    if numareas==5
                        n_m1_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,area5,:))))/length_diss1; 
                    elseif numareas==6
                        n_m1_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,area5,:))))/length_diss1; 
                        n_m1_areas(:,6,sub) = sum(sum(isnan(scan(1:length_diss1,area6,:))))/length_diss1; 
                    end
                end
            else
                for sub=1:2
                    if sub==1
                        scan=deoxy3D(2).sub1;
                    else
                        scan=deoxy3D(2).sub2;
                    end
                    n_m2_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss2,area1,:))))/length_diss2;
                    n_m2_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss2,area2,:))))/length_diss2;
                    n_m2_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss2,area3,:))))/length_diss2;
                    n_m2_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss2,area4,:))))/length_diss2; 
                    if numareas==5
                        n_m2_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss2,area5,:))))/length_diss2; 
                    elseif numareas==6
                        n_m2_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss2,area5,:))))/length_diss2; 
                        n_m2_areas(:,6,sub) = sum(sum(isnan(scan(1:length_diss2,area6,:))))/length_diss2; 
                    end
                end
            end
        end
    end

    %Creates a mask for missing channels for both subjects within each dyad
    for convo=1:2
        if convo==1
            convoMask=n_m1_areas;
        else
            convoMask=n_m2_areas;
        end
        
        for dyad=1:numdyads
            for area=1:numareas 
                if numareas==5 %for all areas
                    if area==1 || area==2
                        if convoMask(dyad,area,1)>2 || convoMask(dyad,area,2)>2
                            nmask(dyad,area,convo)=1;
                        else
                            nmask(dyad,area,convo)=0;
                        end                
                    elseif area==4 || area==5
                        if convoMask(dyad,area,1)>3 || convoMask(dyad,area,2)>3
                            nmask(dyad,area,convo)=1;
                        else
                            nmask(dyad,area,convo)=0;
                        end
                    else
                        if convoMask(dyad,area,1)>1 || convoMask(dyad,area,2)>1
                            nmask(dyad,area,convo)=1;
                        else
                            nmask(dyad,area,convo)=0;
                        end
                    end
                else   % for lateralization or subareas               
                    if area==4 || area==5
                        if convoMask(dyad,area,1)>1 || convoMask(dyad,area,2)>1
                            nmask(dyad,area,convo)=1;
                        else
                            nmask(dyad,area,convo)=0;
                        end
                    else
                        if convoMask(dyad,area,1)>=1 || convoMask(dyad,area,2)>=1
                            nmask(dyad,area,convo)=1;
                        else
                            nmask(dyad,area,convo)=0;
                        end
                    end
                end
            end
        end
    end
    nmask=logical(nmask); %Make it a logical for later masking
    
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

    r_d1=r_d1_areas;
    r_d1(nmask(:,:,1))=NaN;
    Sig_r_d1=r_d1;
    Sig_r_d1(~mask_d1)=NaN;
    Sig_r_d1(:,:,2)=r_d1;
    Sig_r_d1(:,:,3)=r_d1_areas;

    r_d2=r_d2_areas;
    r_d2(nmask(:,:,1))=NaN;
    Sig_r_d2=r_d2;
    Sig_r_d2(~mask_d2)=NaN;
    Sig_r_d2(:,:,2)=r_d2;
    Sig_r_d2(:,:,3)=r_d2_areas;

    if writeXL
        %z-scored for later comparison
        Sig_r_d1=atanh(Sig_r_d1);
        Sig_r_d2=atanh(Sig_r_d2);

        load(strcat('SS_NIRS',filesep,'SS_dyads.mat'))
        dyads=table2array(dyads);
        
        if oxyOnly %file naming for type of oxy
            typeOxy='oxy_';
        else
            typeOxy='deoxy_';
        end
        
        if numareas==4 %variable names and file name for # of areas
            areas=["Dyads","R-mPFC","R-lPFC","R-sms","R-tpj"];
            aName='RLat_';
        elseif numareas==5
            areas=["Dyads","mPFC","lPFC","pmc","sms","tpj"];
            aName='OG_';
        else
            areas=["Dyads","vmPFC","dmPFC","vlPFC","sms","atpj","pTPJ"];
            aName='VD_';
        end
        
        Sig_r_D1=array2table([dyads,Sig_r_d1(:,:,1)],'VariableNames',areas);
        r_d1_lost=array2table([dyads,Sig_r_d1(:,:,2)],'VariableNames',areas);
        r_values_d1=array2table([dyads,Sig_r_d1(:,:,3)],'VariableNames',areas);

        Sig_r_D2=array2table([dyads,Sig_r_d2(:,:,1)],'VariableNames',areas);
        r_d2_lost=array2table([dyads,Sig_r_d2(:,:,2)],'VariableNames',areas);
        r_values_d2=array2table([dyads,Sig_r_d2(:,:,3)],'VariableNames',areas);
        
        if FDR
            xlName=strcat(preprocess_dir,filesep,aName,typeOxy,'D1_mask_',num2str(cutoff),'.xls');
            writetable(r_values_d1,xlName,'sheet','r_vals')
            writetable(r_d1_lost,xlName,'sheet','r_vals_lostchs')
            writetable(Sig_r_D1,xlName,'sheet','r_vals_mask')

            xlName=strcat(preprocess_dir,filesep,aName,typeOxy,'D2_mask_',num2str(cutoff),'.xls');
            writetable(r_values_d2,xlName,'sheet','r_vals')
            writetable(r_d2_lost,xlName,'sheet','r_vals_lostchs')
            writetable(Sig_r_D2,xlName,'sheet','r_vals_mask')
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
    clearvars -except preprocess_dir numchans numdyads image
    
%% Imagine NIRS results
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