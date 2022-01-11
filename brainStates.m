%% Set the Directory
% You can set the directory w/ line 3 or use the get directory in line 4
% preprocess_dir= ''; %Make sure to specify correct location
preprocess_dir=uigetdir('','Choose Data Directory');

%% Analyses to run
% Set to zero if you do not want to perform the task
compile=1;      % If the data has yet to be compiled into all dyad matices
trim=1;         % Trim data to shortest scan
oxyOnly=1;      % 0=z_deoxy; 1=z_oxy
SS=0;           % 0=run conflict study; 1=run social support study
kMeans=0;
%% Properties 
if SS
    % These can be changed based on the study and what FDR cutoff is prefered
    dataprefix='SS';
    ch_reject=3;  % 1=no channel rejection, 2=reject noisy, 3=reject noisy & uncertain
    length_rest=434;
    length_loc1=1513;
    length_loc2=541;
    length_diss1=1853; % Can be changed to the shortest conversation
    length_diss2=1810; % Can be changed to the shortest conversation
    numdyads=52; % Number of dyads
    numchans=42; % Number of channels
    sampRate=3.906250;
else
    % These can be changed based on the study and what FDR cutoff is prefered
    dataprefix='0';
    ch_reject=3;  % 1=no channel rejection, 2=reject noisy, 3=reject noisy & uncertain
    length_rest1=1547;
    length_rest2=1522;
    length_rest3=1549;
    length_affil=2442; % Can be changed to the shortest conversation
    length_con=2468; % Can be changed to the shortest conversation
    numdyads=54; % Number of dyads
    numchans=42; % Number of channels
    sampRate=5;
end

%%
if compile
    numScans=5;
    [deoxy3D,oxy3D] = compileNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans);

    if SS
        %Save a .mat file to the preprocessing directory 
        save(strcat(preprocess_dir,filesep,'SS_compiled'),'deoxy3D','oxy3D');
    else
        %Save a .mat file to the preprocessing directory 
        save(strcat(preprocess_dir,filesep,'CF_compiled'),'deoxy3D','oxy3D');
    end
end

%%
if trim
    if oxyOnly
        if SS
            %Trim to the shortest scan
            oxy3D(1).sub1=oxy3D(1).sub1(1:length_diss1,:,:);
            oxy3D(1).sub2=oxy3D(1).sub2(1:length_diss1,:,:);

            oxy3D(2).sub1=oxy3D(2).sub1(1:length_diss2,:,:);
            oxy3D(2).sub2=oxy3D(2).sub2(1:length_diss2,:,:);

            oxy3D(3).sub1=oxy3D(3).sub1(1:length_loc1,:,:);
            oxy3D(3).sub2=oxy3D(3).sub2(1:length_loc1,:,:);
            
            oxy3D(4).sub1=oxy3D(4).sub1(1:length_loc2,:,:);
            oxy3D(4).sub2=oxy3D(4).sub2(1:length_loc2,:,:);
            
            oxy3D(5).sub1=oxy3D(5).sub1(1:length_rest,:,:);
            oxy3D(5).sub2=oxy3D(5).sub2(1:length_rest,:,:);

        else
                %Trim to the shortest scan
            oxy3D(1).sub1=oxy3D(1).sub1(1:length_affil,:,:);
            oxy3D(1).sub2=oxy3D(1).sub2(1:length_affil,:,:);
            
            oxy3D(2).sub1=oxy3D(2).sub1(1:length_con,:,:);
            oxy3D(2).sub2=oxy3D(2).sub2(1:length_con,:,:);
            
            oxy3D(3).sub1=oxy3D(3).sub1(1:length_rest1,:,:);
            oxy3D(3).sub2=oxy3D(3).sub2(1:length_rest1,:,:);
            
            oxy3D(4).sub1=oxy3D(4).sub1(1:length_rest2,:,:);
            oxy3D(4).sub2=oxy3D(4).sub2(1:length_rest2,:,:);
            
            oxy3D(5).sub1=oxy3D(5).sub1(1:length_rest3,:,:);
            oxy3D(5).sub2=oxy3D(5).sub2(1:length_rest3,:,:);

            oxy3D(5).sub1(:,:,31)=NaN; %Missing scans
            oxy3D(5).sub2(:,:,31)=NaN;
        end
    else
        if SS
            %Trim to the shortest scan
            deoxy3D(1).sub1=deoxy3D(1).sub1(1:length_diss1,:,:);
            deoxy3D(1).sub2=deoxy3D(1).sub2(1:length_diss1,:,:);

            deoxy3D(2).sub1=deoxy3D(2).sub1(1:length_diss2,:,:);
            deoxy3D(2).sub2=deoxy3D(2).sub2(1:length_diss2,:,:);

            deoxy3D(3).sub1=deoxy3D(3).sub1(1:length_loc1,:,:);
            deoxy3D(3).sub2=deoxy3D(3).sub2(1:length_loc1,:,:);
            
            deoxy3D(4).sub1=deoxy3D(4).sub1(1:length_loc2,:,:);
            deoxy3D(4).sub2=deoxy3D(4).sub2(1:length_loc2,:,:);
            
            deoxy3D(5).sub1=deoxy3D(5).sub1(1:length_rest,:,:);
            deoxy3D(5).sub2=deoxy3D(5).sub2(1:length_rest,:,:);
        else
                %Trim to the shortest scan
            deoxy3D(1).sub1=deoxy3D(1).sub1(1:length_affil,:,:);
            deoxy3D(1).sub2=deoxy3D(1).sub2(1:length_affil,:,:);
            
            deoxy3D(2).sub1=deoxy3D(2).sub1(1:length_con,:,:);
            deoxy3D(2).sub2=deoxy3D(2).sub2(1:length_con,:,:);
            
            deoxy3D(3).sub1=deoxy3D(3).sub1(1:length_rest1,:,:);
            deoxy3D(3).sub2=deoxy3D(3).sub2(1:length_rest1,:,:);
            
            deoxy3D(4).sub1=deoxy3D(4).sub1(1:length_rest2,:,:);
            deoxy3D(4).sub2=deoxy3D(4).sub2(1:length_rest2,:,:);
            
            deoxy3D(5).sub1=deoxy3D(5).sub1(1:length_rest3,:,:);
            deoxy3D(5).sub2=deoxy3D(5).sub2(1:length_rest3,:,:);

            deoxy3D(5).sub1(:,:,31)=NaN; %Missing scans
            deoxy3D(5).sub2(:,:,31)=NaN;
        end
    end
end

%%
if oxyOnly
    if SS
        z_oxy_diss1=oxy3D(1).sub1;
        z_oxy_diss1(:,:,53:104)=oxy3D(1).sub2;

        z_oxy_diss2=oxy3D(2).sub1;
        z_oxy_diss2(:,:,53:104)=oxy3D(2).sub2;

        z_oxy_loc1=oxy3D(3).sub1;
        z_oxy_loc1(:,:,53:104)=oxy3D(3).sub2;

        z_oxy_loc2=oxy3D(4).sub1;
        z_oxy_loc2(:,:,53:104)=oxy3D(4).sub2;

        z_oxy_rest=oxy3D(5).sub1;
        z_oxy_rest(:,:,53:104)=oxy3D(5).sub2;
    else
        z_oxy_affil=oxy3D(1).sub1;
        z_oxy_affil(:,:,55:108)=oxy3D(1).sub2;

        z_oxy_con=oxy3D(2).sub1;
        z_oxy_con(:,:,55:108)=oxy3D(2).sub2;

        z_oxy_rest1=oxy3D(3).sub1;
        z_oxy_rest1(:,:,55:108)=oxy3D(3).sub2;

        z_oxy_rest2=oxy3D(4).sub1;
        z_oxy_rest2(:,:,55:108)=oxy3D(4).sub1;

        z_oxy_rest3=oxy3D(5).sub1;
        z_oxy_rest3(:,:,55:108)=oxy3D(5).sub2;
        
        z_oxy(1).name='affil';
        z_oxy(1).scan=z_oxy_affil;
        z_oxy(2).name='con';
        z_oxy(2).scan=z_oxy_con;
        z_oxy(3).name='rest1';
        z_oxy(3).scan=z_oxy_rest1;
        z_oxy(4).name='rest2';
        z_oxy(4).scan=z_oxy_rest2;
        z_oxy(5).name='rest3';
        z_oxy(5).scan=z_oxy_rest3;
    end
else
    if SS
        z_deoxy_diss1=deoxy3D(1).sub1;
        z_deoxy_diss1(:,:,53:104)=deoxy3D(1).sub2;

        z_deoxy_diss2=deoxy3D(2).sub1;
        z_deoxy_diss2(:,:,53:104)=deoxy3D(2).sub2;

        z_deoxy_loc1=deoxy3D(3).sub1;
        z_deoxy_loc1(:,:,53:104)=deoxy3D(3).sub2;

        z_deoxy_loc2=deoxy3D(4).sub1;
        z_deoxy_loc2(:,:,53:104)=deoxy3D(4).sub2;

        z_deoxy_rest=deoxy3D(5).sub1;
        z_deoxy_rest(:,:,53:104)=deoxy3D(5).sub2;  
    else
        z_deoxy_affil=deoxy3D(1).sub1;
        z_deoxy_affil(:,:,55:108)=deoxy3D(1).sub2;

        z_deoxy_con=deoxy3D(2).sub1;
        z_deoxy_con(:,:,55:108)=deoxy3D(2).sub2;

        z_deoxy_rest1=deoxy3D(3).sub1;
        z_deoxy_rest1(:,:,55:108)=deoxy3D(3).sub2;

        z_deoxy_rest2=deoxy3D(4).sub1;
        z_deoxy_rest2(:,:,55:108)=deoxy3D(4).sub1;

        z_deoxy_rest3=deoxy3D(5).sub1;
        z_deoxy_rest3(:,:,55:108)=deoxy3D(5).sub2;
    end
end

%%
if kMeans
    %Reduce the dimensionality for k-means clustering
    z_oxyAffil=z_oxy(1).scan(:,:,sub);
    z_oxyCon=z_oxy(2).scan(:,:,sub);
    z_oxyRest1=z_oxy(3).scan(:,:,sub);
    z_oxyRest2=z_oxy(4).scan(:,:,sub);
    z_oxyRest3=z_oxy(5).scan(:,:,sub);
    for sub=2:(numdyads*2)
        z_oxyAffil=[z_oxyAffil;z_oxy(1).scan(:,:,sub)];
        z_oxyCon=[z_oxyCon;z_oxy(2).scan(:,:,sub)];
        z_oxyRest1=[z_oxyRest1;z_oxy(3).scan(:,:,sub)];
        z_oxyRest2=[z_oxyRest2;z_oxy(4).scan(:,:,sub)];
        z_oxyRest3=[z_oxyRest3;z_oxy(5).scan(:,:,sub)];
    end

%     for ch=1:numchans
%         for tp=1:length(z_oxyAffil)
%             if isnan(z_oxyAffil(tp,ch))
%                 z_oxyAffil(tp,ch)=0;
%             end
%         end
%     end
            
    figure;
    plot(z_oxyAffil(:,1),z_oxyAffil(:,8),'k*','MarkerSize',2);

    %LEFT OFF HERE
    rng(1); % For reproducibility
    [idx,C] = kmeans(z_oxyAffil,5);

    x1 = min(z_oxyAffil(:,1)):0.01:max(z_oxyAffil(:,1));
    x2 = min(z_oxyAffil(:,2)):0.01:max(z_oxyAffil(:,2));
    [x1G,x2G] = meshgrid(x1,x2);
    XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot

    idx2Region = kmeans(XGrid,5,'MaxIter',1,'Start',C);

    figure;
    gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
        [0,0.75,0.75;0.75,0,0.75;0.75,0.75,0],'..');
    hold on;
    plot(X(:,1),X(:,2),'k*','MarkerSize',5);
    title 'Fisher''s Iris Data';
    xlabel 'Petal Lengths (cm)';
    ylabel 'Petal Widths (cm)'; 
    legend('Region 1','Region 2','Region 3','Data','Location','SouthEast');
    hold off;
end
    