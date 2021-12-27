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
compile=0;      % If the data has yet to be compiled into all dyad matices
oxyOnly=1;      % 0=z_totaloxy; 1=z_oxy
chCorr=0;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=1;     % Dyadic correlation over entire conversation per area of the brain
FDR=1;          % False discovery rate correction
writeXL=1;      % If you want to write the data to an excel sheet(s)
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
numareas=5; % Number of areas of the brain. Must match the combined ch's in the 'areaCorr' section

%% Compile subject data %%
% This will compile the two conversations of interest. Creates a 3D
% matrix for type of oxy (deoxy,oxy,totaloxy), discussion (1 and 2),  
% and subject (1 and 2) for twelve (3x2x2=12) 3D matrices (time,channels,dyad)
if compile
    [z_deoxy1_1,z_oxy1_1,z_totaloxy1_1,z_deoxy1_2,z_oxy1_2,z_totaloxy1_2,...
        z_deoxy2_1,z_oxy2_1,z_totaloxy2_1,z_deoxy2_2,z_oxy2_2,...
        z_totaloxy2_2] = compileNIRSdata(preprocess_dir,dataprefix,ch_reject);

    z_deoxy1_diss1=z_deoxy1_1;
    z_deoxy2_diss1=z_deoxy1_2;
    z_deoxy1_diss2=z_deoxy2_1;
    z_deoxy2_diss2=z_deoxy2_2;
    z_oxy1_diss1=z_oxy1_1;
    z_oxy2_diss1=z_oxy1_2;
    z_oxy1_diss2=z_oxy2_1;
    z_oxy2_diss2=z_oxy2_2;
    z_totaloxy1_diss1=z_totaloxy1_1;
    z_totaloxy2_diss1=z_totaloxy1_2;
    z_totaloxy1_diss2=z_totaloxy2_1;
    z_totaloxy2_diss2=z_totaloxy2_2;
    %Save a .mat file to the preprocessing directory 
    save(strcat(preprocess_dir,filesep,'SS_compiled'),'z_deoxy1_diss1','z_deoxy2_diss1',...
    'z_oxy1_diss1','z_oxy2_diss1','z_totaloxy1_diss1','z_totaloxy2_diss1',...
    'z_deoxy1_diss2','z_deoxy2_diss2','z_oxy1_diss2','z_oxy2_diss2',...
    'z_totaloxy1_diss2','z_totaloxy2_diss2');

    clearvars -except preprocess_dir numdyads numchans length_diss1 length_diss2 oxyOnly chCorr areaCorr cutoff FDR image writeXL
end

%% Compute basic channel correlations %%
% Computes the matched channel correlations for each dyad and conversation. 
if chCorr
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'z_oxy1_diss1',...
            'z_oxy2_diss1','z_oxy1_diss2','z_oxy2_diss2');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = z_oxy1_diss1(1:length_diss1,channel,dyad);
                b = z_oxy2_diss1(1:length_diss1,channel,dyad);
                [r_values_diss1(dyad,channel),p_values_diss1(dyad,channel)] = corr(a,b);

                c = z_oxy1_diss2(1:length_diss2,channel,dyad);
                d = z_oxy2_diss2(1:length_diss2,channel,dyad);
                [r_values_diss2(dyad,channel),p_values_diss2(dyad,channel)] = corr(c,d);
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'z_totaloxy1_diss1',...
            'z_totaloxy2_diss1','z_totaloxy1_diss2','z_totaloxy2_diss2');
        
        for dyad=1:numdyads
            for channel=1:numchans
                a = z_totaloxy1_diss1(:,channel,dyad);
                b = z_totaloxy2_diss1(:,channel,dyad);
                [r_values_diss1(dyad,channel),p_values_diss1(dyad,channel)] = corr(a,b);

                c = z_totaloxy1_diss2(:,channel,dyad);
                d = z_totaloxy2_diss2(:,channel,dyad);
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

    clearvars -except preprocess_dir numdyads numchans length_diss1 length_diss2 oxyOnly chCorr areaCorr cutoff FDR image writeXL
end

%% Mean Timecourse Synchrony per ROI per Dyad %%
% Averages activiation in area specific channels then computes dyadic 
% correlation of those areas of the brain 
if areaCorr   
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'z_oxy1_diss1',...
            'z_oxy2_diss1','z_oxy1_diss2','z_oxy2_diss2');
        
        %Get the mean for each area of interest. 
        %Number of channels:mpfc[8 chs],lpfc[6 chs],pmc[4 chs],sms[10chs],tpj[12 chs]
        z1_diss1_areas(:,:,1) = nanmean(z_oxy1_diss1(1:length_diss1,[1:3,5,8,10:12],:),2);  %medial prefrontal cortex
        z1_diss1_areas(:,:,2) = nanmean(z_oxy1_diss1(1:length_diss1,[4,6,7,9,30:31],:),2);  %lateral prefrontal cortex
        z1_diss1_areas(:,:,3) = nanmean(z_oxy1_diss1(1:length_diss1,[13,24:25,32],:),2);     %premotor/motor cortex
        z1_diss1_areas(:,:,4) = nanmean(z_oxy1_diss1(1:length_diss1,[14:16,26:29,33:35],:),2);%somatosensory cortex
        z1_diss1_areas(:,:,5) = nanmean(z_oxy1_diss1(1:length_diss1,[17:21,23,36:40,42],:),2);%temporoparietal junction
        
        z2_diss1_areas(:,:,1) = nanmean(z_oxy2_diss1(1:length_diss1,[1:3,5,8,10:12],:),2); 
        z2_diss1_areas(:,:,2) = nanmean(z_oxy2_diss1(1:length_diss1,[4,6,7,9,30:31],:),2);
        z2_diss1_areas(:,:,3) = nanmean(z_oxy2_diss1(1:length_diss1,[13,24:25,32],:),2);
        z2_diss1_areas(:,:,4) = nanmean(z_oxy2_diss1(1:length_diss1,[14:16,26:29,33:35],:),2);
        z2_diss1_areas(:,:,5) = nanmean(z_oxy2_diss1(1:length_diss1,[17:21,23,36:40,42],:),2);
        
        z1_diss2_areas(:,:,1) = nanmean(z_oxy1_diss2(1:length_diss2,[1:3,5,8,10:12],:),2);
        z1_diss2_areas(:,:,2) = nanmean(z_oxy1_diss2(1:length_diss2,[4,6,7,9,30:31],:),2);
        z1_diss2_areas(:,:,3) = nanmean(z_oxy1_diss2(1:length_diss2,[13,24:25,32],:),2);
        z1_diss2_areas(:,:,4) = nanmean(z_oxy1_diss2(1:length_diss2,[14:16,26:29,33:35],:),2);
        z1_diss2_areas(:,:,5) = nanmean(z_oxy1_diss2(1:length_diss2,[17:21,23,36:40,42],:),2);
        
        z2_diss2_areas(:,:,1) = nanmean(z_oxy2_diss2(1:length_diss2,[1:3,5,8,10:12],:),2); 
        z2_diss2_areas(:,:,2) = nanmean(z_oxy2_diss2(1:length_diss2,[4,6,7,9,30:31],:),2);
        z2_diss2_areas(:,:,3) = nanmean(z_oxy2_diss2(1:length_diss2,[13,24:25,32],:),2);
        z2_diss2_areas(:,:,4) = nanmean(z_oxy2_diss2(1:length_diss2,[14:16,26:29,33:35],:),2);
        z2_diss2_areas(:,:,5) = nanmean(z_oxy2_diss2(1:length_diss2,[17:21,23,36:40,42],:),2);
        
        %Number of missing channels per subject, per area of interest
        for dc=1:2
            if dc==1
                for sub=1:2
                    if sub==1
                        scan=z_oxy1_diss1;
                    else
                        scan=z_oxy2_diss1;
                    end
                    n_m1_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss1,[1:3,5,8,10:12],:))))/length_diss1;
                    n_m1_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss1,[4,6,7,9,30:31],:))))/length_diss1;
                    n_m1_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss1,[13,24:25,32],:))))/length_diss1;
                    n_m1_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss1,[14:16,26:29,33:35],:))))/length_diss1;        
                    n_m1_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,[17:21,36:40],:))))/length_diss1;
                end
            else
                for sub=1:2
                    if sub==1
                        scan=z_oxy1_diss2;
                    else
                        scan=z_oxy2_diss2;
                    end
                    n_m2_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss1,[1:3,5,8,10:12],:))))/length_diss2;
                    n_m2_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss1,[4,6,7,9,30:31],:))))/length_diss2;
                    n_m2_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss1,[13,24:25,32],:))))/length_diss2;
                    n_m2_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss1,[14:16,26:29,33:35],:))))/length_diss2;        
                    n_m2_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,[17:21,36:40],:))))/length_diss2;
                end
            end
        end
      
    else
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'z_totaloxy1_diss1',...
            'z_totaloxy2_diss1','z_totaloxy1_diss2','z_totaloxy2_diss2');
         
        %Get the mean for each area of interest. 
        %Number of channels:mpfc[8 chs],lpfc[6 chs],pmc[4 chs],sms[10chs],tpj[12 chs]
        z1_diss1_areas(:,:,1) = nanmean(z_totaloxy1_diss1(1:length_diss1,[1:3,5,8,10:12],:),2); %mpfc
        z1_diss1_areas(:,:,2) = nanmean(z_totaloxy1_diss1(1:length_diss1,[4,6,7,9,30:31],:),2); %lpfc
        z1_diss1_areas(:,:,3) = nanmean(z_totaloxy1_diss1(1:length_diss1,[13,24:25,32],:),2);  %pmc
        z1_diss1_areas(:,:,4) = nanmean(z_totaloxy1_diss1(1:length_diss1,[14:16,26:29,33:35],:),2); %sms        
        z1_diss1_areas(:,:,5) = nanmean(z_totaloxy1_diss1(1:length_diss1,[17:21,36:40],:),2); %tpj
        
        z2_diss1_areas(:,:,1) = nanmean(z_totaloxy2_diss1(1:length_diss1,[1:3,5,8,10:12],:),2); 
        z2_diss1_areas(:,:,2) = nanmean(z_totaloxy2_diss1(1:length_diss1,[4,6,7,9,30:31],:),2);
        z2_diss1_areas(:,:,3) = nanmean(z_totaloxy2_diss1(1:length_diss1,[13,24:25,32],:),2);
        z2_diss1_areas(:,:,4) = nanmean(z_totaloxy2_diss1(1:length_diss1,[14:16,26:29,33:35],:),2);        
        z2_diss1_areas(:,:,5) = nanmean(z_totaloxy2_diss1(1:length_diss1,[17:21,36:40],:),2);
        
        z1_diss2_areas(:,:,1) = nanmean(z_totaloxy1_diss2(1:length_diss2,[1:3,5,8,10:12],:),2);
        z1_diss2_areas(:,:,2) = nanmean(z_totaloxy1_diss2(1:length_diss2,[4,6,7,9,30:31],:),2);
        z1_diss2_areas(:,:,3) = nanmean(z_totaloxy1_diss2(1:length_diss2,[13,24:25,32],:),2);
        z1_diss2_areas(:,:,4) = nanmean(z_totaloxy1_diss2(1:length_diss2,[14:16,26:29,33:35],:),2);        
        z1_diss2_areas(:,:,5) = nanmean(z_totaloxy1_diss2(1:length_diss2,[17:21,36:40],:),2);
        
        z2_diss2_areas(:,:,1) = nanmean(z_totaloxy2_diss2(1:length_diss2,[1:3,5,8,10:12],:),2); 
        z2_diss2_areas(:,:,2) = nanmean(z_totaloxy2_diss2(1:length_diss2,[4,6,7,9,30:31],:),2);
        z2_diss2_areas(:,:,3) = nanmean(z_totaloxy2_diss2(1:length_diss2,[13,24:25,32],:),2);
        z2_diss2_areas(:,:,4) = nanmean(z_totaloxy2_diss2(1:length_diss2,[14:16,26:29,33:35],:),2);        
        z2_diss2_areas(:,:,5) = nanmean(z_totaloxy2_diss2(1:length_diss2,[17:21,36:40],:),2);
    
    %Number of missing channels per subject, per area of interest
        for dc=1:2
            if dc==1
                for sub=1:2
                    if sub==1
                        scan=z_totaloxy1_diss1;
                    else
                        scan=z_totaloxy2_diss1;
                    end
                    n_m1_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss1,[1:3,5,8,10:12],:))))/length_diss1;
                    n_m1_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss1,[4,6,7,9,30:31],:))))/length_diss1;
                    n_m1_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss1,[13,24:25,32],:))))/length_diss1;
                    n_m1_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss1,[14:16,26:29,33:35],:))))/length_diss1;        
                    n_m1_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,[17:21,36:40],:))))/length_diss1;
                end
            else
                for sub=1:2
                    if sub==1
                        scan=z_totaloxy1_diss2;
                    else
                        scan=z_totaloxy2_diss2;
                    end
                    n_m2_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss2,[1:3,5,8,10:12],:))))/length_diss2;
                    n_m2_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss2,[4,6,7,9,30:31],:))))/length_diss2;
                    n_m2_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss2,[13,24:25,32],:))))/length_diss2;
                    n_m2_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss2,[14:16,26:29,33:35],:))))/length_diss2;        
                    n_m2_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss2,[17:21,36:40],:))))/length_diss2;
                end
            end
        end
    end

    %Creates a mask for missing channels based on both subjects within each
    %dyad. Anything over 1/4 missing channels mask as 0. Number of channels:
    %mpfc[8/4=2chs],lpfc[6/4=2chs],pmc[4/4=1chs],sms[10/4=3chs],tpj[12/4=3chs]
    for convo=1:2
        if convo==1
            convoMask=n_m1_areas;
        else
            convoMask=n_m2_areas;
        end
        
        for dyad=1:numdyads
            for area=1:numareas 
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
            end
        end
    end
    nmask=logical(nmask); %Make it a logical for later masking
    
    % Dyadic correlation per Conversation
    for dyad=1:numdyads
        for area=1:numareas
            a = z1_diss1_areas(:,dyad,area);
            b = z2_diss1_areas(:,dyad,area);
            [r_d1_areas(dyad,area),p_d1_areas(dyad,area)] = corr(a,b);

            c = z1_diss2_areas(:,dyad);
            d = z2_diss2_areas(:,dyad);
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

    save(strcat(preprocess_dir,filesep,'SS_areas'),'Sig_r_d1','Sig_r_d2')

    if writeXL
        %z-scored for later comparison
        Sig_r_d1=atanh(Sig_r_d1);
        Sig_r_d2=atanh(Sig_r_d2);

        load(strcat('SS_NIRS',filesep,'SS_dyads.mat'))
        dyads=table2array(dyads);
        areas=["Dyads","mPFC","lPFC","pmc","sms","tpj"];
        
        Sig_r_D1=array2table([dyads,Sig_r_d1(:,:,1)],'VariableNames',areas);
        r_d1_lost=array2table([dyads,Sig_r_d1(:,:,2)],'VariableNames',areas);
        r_values_d1=array2table([dyads,Sig_r_d1(:,:,3)],'VariableNames',areas);

        Sig_r_D2=array2table([dyads,Sig_r_d2(:,:,1)],'VariableNames',areas);
        r_d2_lost=array2table([dyads,Sig_r_d2(:,:,2)],'VariableNames',areas);
        r_values_d2=array2table([dyads,Sig_r_d2(:,:,3)],'VariableNames',areas);

        if FDR
            xlName=strcat(preprocess_dir,filesep,'D1_Sig_mask_',num2str(cutoff),'.xls');
            writetable(r_values_d1,xlName,'sheet','r_vals')
            writetable(r_d1_lost,xlName,'sheet','r_vals_lostchs')
            writetable(Sig_r_D1,xlName,'sheet','r_vals_mask')

            xlName=strcat(preprocess_dir,filesep,'D2_Sig_mask_',num2str(cutoff),'.xls');
            writetable(r_values_d2,xlName,'sheet','r_vals')
            writetable(r_d2_lost,xlName,'sheet','r_vals_lostchs')
            writetable(Sig_r_D2,xlName,'sheet','r_vals_mask')
        else
            xlName=strcat('SS_areas.xls');
            writetable(r_d1_areas,xlName,'sheet','r_d1_areas')
            writetable(p_d1_areas,xlName,'sheet','p_d1_areas')
            writetable(r_d2_areas,xlName,'sheet','r_d2_areas')
            writetable(p_d2_areas,xlName,'sheet','p_d2_areas')
        end
    end
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