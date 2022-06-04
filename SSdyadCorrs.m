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
numareas=4; % 4=mPFC, lPFC, TPJ, VIS; 5=mPFC, lPFC, PMC, SMS & TPJ; 6=vmPFC, dmPFC, vlPFC, SMS, TPJ
zdim=1; % 0=compile non-zscored; 1=compile z-scored

%% Compile subject data %%
% Compiles the two conversations of interest. Creates a 3D
% matrix for type of oxy (deoxy & oxy), discussion (1 and 2),  
% and subject (1 and 2) for eight (2x2x2=12) 3D matrices (time,channels,dyad)
if compile
    numScans=5;
    [deoxy3D,oxy3D]= compiledyadicNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans,zdim);

    %Save a .mat file to the preprocessing directory 
    save(strcat(preprocess_dir,filesep,'SS_compiled'),'deoxy3D','oxy3D');
    
    clearvars -except preprocess_dir numdyads numchans numareas length_diss1 length_diss2 oxyOnly chCorr areaCorr cutoff FDR writeXL
end

%% Compute all channel correlations %%
% Computes all channel correlations for each dyad and conversation. 
if chCorr
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'oxy3D');
        r_values_diss1=nan(numchans,numchans,numdyads);
        r_values_diss2=nan(numchans,numchans,numdyads);
        p_values_diss1=nan(numchans,numchans,numdyads);
        p_values_diss2=nan(numchans,numchans,numdyads);
        
        for dyad=1:numdyads
            for chan1=1:numchans
                for chan2=1:numchans
                    a = oxy3D(1).sub1(1:length_diss1,chan1,dyad);
                    b = oxy3D(1).sub2(1:length_diss1,chan2,dyad);
                    [r_values_diss1(chan1,chan2,dyad),p_values_diss1(chan1,chan2,dyad)] = corr(a,b);
    
                    c = oxy3D(2).sub1(1:length_diss2,chan1,dyad);
                    d = oxy3D(2).sub2(1:length_diss2,chan2,dyad);
                    [r_values_diss2(chan1,chan2,dyad),p_values_diss2(chan1,chan2,dyad)] = corr(c,d);
                end
            end
        end
    else
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'deoxy3D');
        r_values_diss1=nan(numchans,numchans,numdyads);
        r_values_diss2=nan(numchans,numchans,numdyads);
        p_values_diss1=nan(numchans,numchans,numdyads);
        p_values_diss2=nan(numchans,numchans,numdyads);

        for dyad=1:numdyads
            for chan1=1:numchans
                for chan2=1:numchans
                    a = oxy3D(1).sub1(1:length_diss1,chan1,dyad);
                    b = oxy3D(1).sub2(1:length_diss1,chan2,dyad);
                    [r_values_diss1(chan1,chan2,dyad),p_values_diss1(chan1,chan2,dyad)] = corr(a,b);
    
                    c = oxy3D(2).sub1(1:length_diss2,chan1,dyad);
                    d = oxy3D(2).sub2(1:length_diss2,chan2,dyad);
                    [r_values_diss2(chan1,chan2,dyad),p_values_diss2(chan1,chan2,dyad)] = corr(c,d);
                end
            end
        end
    end

    %for paired channels corrs only
    for dy=1:numdyads
        for ch=1:numchans
            r_paired_diss1(ch,dy)=r_values_diss1(ch,ch,dy);
            r_paired_diss2(ch,dy)=r_values_diss2(ch,ch,dy);
        end
    end

    %Uses the FDR_bky correction based on the cutoff specified in properties
    if FDR
        r_mask_diss1=r_values_diss1;
        r_mask_diss2=r_values_diss2;
        for dy=1:numdyads
            [mask_diss1(:,:,dy), ~]=fdr_bky(p_values_diss1(:,:,dy),cutoff,'yes'); 
            [mask_diss2(:,:,dy), ~]=fdr_bky(p_values_diss2(:,:,dy),cutoff,'yes');

            for ch=1:numchans
                r_mask_diss1(~mask_diss1(:,ch,dy),ch,dy)=NaN;
                r_mask_diss2(~mask_diss2(:,ch,dy),ch,dy)=NaN;
            end
        end

        save(strcat(preprocess_dir,filesep,'FDR_r_mask'),'FDR_r_mask')
        save(strcat(preprocess_dir,filesep,'SS_FDR_chCorrs'),'r_values_diss1','r_values_diss2',...
        'p_values_diss1','p_values_diss2','r_mask_diss1','r_mask_diss2', 'r_paired_diss1','r_paired_diss2')
    else
        save(strcat(preprocess_dir,filesep,'SS_rVals_dyadChs'),'r_values_diss1','r_values_diss2',...
        'p_values_diss1','p_values_diss2','r_paired_diss1','r_paired_diss2')
    end

    clearvars -except preprocess_dir numdyads numchans numareas length_diss1 length_diss2 oxyOnly chCorr areaCorr cutoff FDR writeXL
end

%% Mean Timecourse Synchrony per ROI per Dyad %%
% Averages activiation in area specific channels then computes dyadic 
% correlation of those areas of the brain 
if areaCorr    
    montageMatch=readtable(strcat(preprocess_dir,filesep,'montageMatch.csv'));

    if numareas==4 
        areas1={[1:3,5,8,10:12];[4,6:7,9,30:31];[17:21,23,36:40,42];[22,41]}; %mPFC, lPFC, SSC, TPJ
        areas2={[1:3,5,8,10];[4,6:7,9,30:31];[16,17:21,35:40];41:42}; %mPFC, lPFC, TPJ, VIS
        areas3={[1:3,5,8];[4,7,30:31];[14:18,33:37];[21:23,40:42]}; %mPFC, lPFC, TPJ,VIS
    elseif numareas==5
        areas={[1:3,5,8,10:12];[4,6,7,9,30:31];[13,24:25,32];[14:16,26:29,33:35];[17:21,23,36:40,42]}; %mPFC, lPFC, PMC, SSC, TPJ
    elseif numareas==6
        areas={1:3;11:12;[4,7,30:31];[14:16,33:35];[17:19,36:38];[21,23,40,42]};%vmPFC, dmPFC, vlPFC, SSC, aTPJ, pTPJ
    end  
    
    if oxyOnly
        load(strcat(preprocess_dir,filesep,'SS_compiled'),'oxy3D');

        [z1_diss1_areas,z1_diss2_areas,z2_diss1_areas,z2_diss2_areas,missN_s1,missN_s2]=areaMeans(oxy3D,...
            length_diss1,length_diss2,numdyads,numareas,montageMatch,areas1,areas2,areas3); %Adjusted area means by cap placement

        %Creates a mask for missing channels for both subjects within each dyad
        [nmask1,nmask2]=maskmissing(montageMatch,missN_s1,missN_s2,numdyads,numareas);
        nmask=nmask1+nmask2;
        nmask=nmask-1;
        nmask=logical(nmask);
    end
    
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
    r_d1=r_d1_areas;
    r_d1(~nmask(:,:,1))=NaN;
    Sig_r_d1=r_d1;
    Sig_r_d1(~mask_d1)=NaN;
    Sig_r_d1(:,:,2)=r_d1;
    Sig_r_d1(:,:,3)=r_d1_areas;

    r_d2=r_d2_areas;
    r_d2(~nmask(:,:,2))=NaN;
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
            areas=["Dyads","mPFC","lPFC","TPJ","VIS"];
            aName='Adj_';
        elseif numareas==5
            areas=["Dyads","mPFC","lPFC","pmc","sms","tpj"];
            aName='OG_';
        else
            areas=["Dyads","vmPFC","dmPFC","vlPFC","sms","atpj","pTPJ"];
            aName='SA_';
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
    clearvars -except preprocess_dir numchans numdyads
