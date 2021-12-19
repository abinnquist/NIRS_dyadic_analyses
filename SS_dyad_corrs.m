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
compile=0;      % If the data has yet to be compiled into all dyad matices
oxyOnly=0;      % 0=z_totaloxy; 1=z_oxy
chCorr=0;       % Dyadic channel correlations over entire conversation (same channels only)
areaCorr=0;     % Dyadic correlation over entire conversation per area of the brain
FDR=0;          % False discovery rate correction
image=0;        % To image a conversation & dyad (will prompt you for which ones)
writeXL=0;      % If you want to write the data to an excel sheet(s)

%% Set the Directory
% You can set the directory w/ line 3 or use the get directory in line 4
% preprocess_dir= ''; %Make sure to specify correct location
preprocess_dir=uigetdir('','Choose Data Directory');

%% Properties
% These can be changed based on the study and what FDR cutoff is prefered
dataprefix='SS';
cutoff=0.1; %P-value to use as false discovery rate cut-off
length_diss1=1853; %Can be changed to the shortest conversation
length_diss2=1810; %Can be changed to the shortest conversation
numdyads=52; 
numchans=42;
numareas=5;

%% Compile subject data %%
% This will compile the two conversations of interest. Creates a 3D
% matrix for type of oxy (deoxy,oxy,totaloxy), discussion (1 and 2),  
% and subject (1 and 2) for twelve (3x2x2=12) 3D matrices (time,channels,dyad)
if compile
[z_deoxy1_1,z_oxy1_1,z_totaloxy1_1,z_deoxy1_2,z_oxy1_2,z_totaloxy1_2,z_deoxy2_1,...
    z_oxy2_1,z_totaloxy2_1,z_deoxy2_2,z_oxy2_2,z_totaloxy2_2] = compileNIRSdata(preprocess_dir,dataprefix);

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
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'z_totaloxy1_diss1_all',...
            'z_totaloxy2_diss1_all','z_totaloxy1_diss2_all','z_totaloxy2_diss2_all');
        
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
        save(strcat(preprocess_dir,filesep,'SS_FDR_chCorrs'),'r_mask_diss1','r_mask_diss2','r_values_diss1','r_values_diss2',...
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
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'z_oxy1_diss1_all',...
            'z_oxy2_diss1_all','z_oxy1_diss2_all','z_oxy2_diss2_all');
        
        %Get the mean for each area of interest. 
        %Number of channels:mpfc[8 chs],lpfc[6 chs],pmc[4 chs],sms[10chs],tpj[12 chs]
        z1_diss1_areas(:,:,1) = nanmean(z_oxy1_diss1_all(1:length_diss1,[1:3,5,8,10:12],:),2);  %medial prefrontal cortex
        z1_diss1_areas(:,:,1) = nanmean(z_oxy1_diss1_all(1:length_diss1,[4,6,7,9,30:31],:),2);  %lateral prefrontal cortex
        z1_diss1_areas(:,:,1) = nanmean(z_oxy1_diss1_all(1:length_diss1,[13,24:25,32],:),2);     %premotor/motor cortex
        z1_diss1_areas(:,:,1) = nanmean(z_oxy1_diss1_all(1:length_diss1,[14:16,26:29,33:35],:),2);%somatosensory cortex
        z1_diss1_areas(:,:,1) = nanmean(z_oxy1_diss1_all(1:length_diss1,[17:21,23,36:40,42],:),2);%temporoparietal junction
        
        z2_diss1_areas(:,:,1) = nanmean(z_oxy2_diss1_all(1:length_diss1,[1:3,5,8,10:12],:),2); 
        z2_diss1_areas(:,:,1) = nanmean(z_oxy2_diss1_all(1:length_diss1,[4,6,7,9,30:31],:),2);
        z2_diss1_areas(:,:,1) = nanmean(z_oxy2_diss1_all(1:length_diss1,[13,24:25,32],:),2);
        z2_diss1_areas(:,:,1) = nanmean(z_oxy2_diss1_all(1:length_diss1,[14:16,26:29,33:35],:),2);
        z2_diss1_areas(:,:,1) = nanmean(z_oxy2_diss1_all(1:length_diss1,[17:21,23,36:40,42],:),2);
        
        z1_diss2_areas(:,:,1) = nanmean(z_oxy1_diss2(1:length_diss2,[1:3,5,8,10:12],:),2);
        z1_diss2_areas(:,:,1) = nanmean(z_oxy1_diss2(1:length_diss2,[4,6,7,9,30:31],:),2);
        z1_diss2_areas(:,:,1) = nanmean(z_oxy1_diss2(1:length_diss2,[13,24:25,32],:),2);
        z1_diss2_areas(:,:,1) = nanmean(z_oxy1_diss2(1:length_diss2,[14:16,26:29,33:35],:),2);
        z1_diss2_areas(:,:,1) = nanmean(z_oxy1_diss2(1:length_diss2,[17:21,23,36:40,42],:),2);
        
        z2_diss2_areas(:,:,1) = nanmean(z_oxy2_diss2(1:length_diss2,[1:3,5,8,10:12],:),2); 
        z2_diss2_areas(:,:,1) = nanmean(z_oxy2_diss2(1:length_diss2,[4,6,7,9,30:31],:),2);
        z2_diss2_areas(:,:,1) = nanmean(z_oxy2_diss2(1:length_diss2,[13,24:25,32],:),2);
        z2_diss2_areas(:,:,1) = nanmean(z_oxy2_diss2(1:length_diss2,[14:16,26:29,33:35],:),2);
        z2_diss2_areas(:,:,1) = nanmean(z_oxy2_diss2(1:length_diss2,[17:21,23,36:40,42],:),2);
        
        %Number of missing channels per subject, per area of interest
        for dc=1:2
            if dc==1
                for sub=1:2
                    if sub==1
                        scan=z_oxy1_diss1_all;
                    else
                        scan=z_oxy2_diss1_all;
                    end
                    n_diss1_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss1,[1:3,5,8,10:12],:))))/length_diss1;
                    n_diss1_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss1,[4,6,7,9,30:31],:))))/length_diss1;
                    n_diss1_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss1,[13,24:25,32],:))))/length_diss1;
                    n_diss1_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss1,[14:16,26:29,33:35],:))))/length_diss1;        
                    n_diss1_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,[17:21,36:40],:))))/length_diss1;
                end
            else
                for sub=1:2
                    if sub==1
                        scan=z_oxy1_diss2;
                    else
                        scan=z_oxy2_diss2;
                    end
                    n_diss2_areas(:,1,sub) = sum(sum(isnan(scan(1:length_diss1,[1:3,5,8,10:12],:))))/length_diss2;
                    n_diss2_areas(:,2,sub) = sum(sum(isnan(scan(1:length_diss1,[4,6,7,9,30:31],:))))/length_diss2;
                    n_diss2_areas(:,3,sub) = sum(sum(isnan(scan(1:length_diss1,[13,24:25,32],:))))/length_diss2;
                    n_diss2_areas(:,4,sub) = sum(sum(isnan(scan(1:length_diss1,[14:16,26:29,33:35],:))))/length_diss2;        
                    n_diss2_areas(:,5,sub) = sum(sum(isnan(scan(1:length_diss1,[17:21,36:40],:))))/length_diss2;
                end
            end
        end
      
    else
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'z_totaloxy1_diss1_all',...
            'z_totaloxy2_diss1_all','z_totaloxy1_diss2_all','z_totaloxy2_diss2_all');
         
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
    
    %Creates mask for missing channels based on both subjects within each
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
    
    if FDR
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

        load(strcat('SS_NIRS',filesep,'SS_dyads.mat'))
        areas=["Dyads","mPFC","lPFC","pmc","sms","tpj"];
        Sig_r_D1=array2table([dyads,Sig_r_d1(:,:,1)],'VariableNames',areas);
        r_d1_lost=array2table([dyads,Sig_r_d1(:,:,2)],'VariableNames',areas);
        r_values_d1=array2table([dyads,Sig_r_d1(:,:,3)],'VariableNames',areas);

        Sig_r_D2=array2table([dyads,Sig_r_d2(:,:,1)],'VariableNames',areas);
        r_d2_lost=array2table([dyads,Sig_r_d2(:,:,2)],'VariableNames',areas);
        r_values_d2=array2table([dyads,Sig_r_d2(:,:,3)],'VariableNames',areas);
    end

    if writeXL
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
            xlName=strcat(preprocess_dir,filesep,'SS_areas.xls');
            writetable(r_d1_areas,xlName,'sheet','r_d1_areas')
            writetable(p_d1_areas,xlName,'sheet','p_d1_areas')
            writetable(r_d2_areas,xlName,'sheet','r_d2_areas')
            writetable(p_d2_areas,xlName,'sheet','p_d2_areas')
        end
    end
end
    clearvars -except preprocess_dir numchans numdyads image

%% Imaging correlated areas
if image
    load(strcat(pwd,filesep,'SS_NIRS',filesep,'SS_mnicoords.mat'));
    load(strcat(preprocess_dir,filesep,'SS_FDR_chCorrs.mat'));

    %Choose the dyad and conversation you wish to image
    convoQ = 'Do you want to image the first (1) disccusion or second (2) Conversation? \n';
    convo = input(convoQ);
    
    if convo==1
        convo2img=r_mask_diss1;
        conversation='diss1';
    else
        convo2img=r_mask_diss2;
        conversation='diss2';
    end
    
    dyadQ = 'What dyad would you like to image (1 to 52)? \n';
    dyad2img = input(dyadQ);

    imgName=strcat('dyad_',num2str(dyad2img),'_',conversation,'.img');

    % Which conversation correlations you want to visualize
    convo_mask=convo2img(dyad2img,:)';

    addpath(genpath(strcat(pwd,filesep,'xjview')))
    addpath(genpath(strcat(pwd,filesep,'spm12')))
    % Make sure to change the name of the file
    nirs2img(imgName, SS_mnicoords, convo_mask, 0,0,0)
end
clear