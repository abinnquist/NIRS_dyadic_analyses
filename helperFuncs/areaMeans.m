function [z1_diss1_areas,z1_diss2_areas,z2_diss1_areas,z2_diss2_areas,missN_s1,missN_s2]=areaMeans(oxy3D,...
    length_diss1,length_diss2,numdyads,numareas,montageMatch,areas1,areas2,areas3)

z1_diss1_areas=nan(length_diss1,numareas,numdyads);
z2_diss1_areas=nan(length_diss1,numareas,numdyads);
z1_diss2_areas=nan(length_diss2,numareas,numdyads);
z2_diss2_areas=nan(length_diss2,numareas,numdyads);

missN_s1=zeros(numdyads,numareas,2);
missN_s2=zeros(numdyads,numareas,2);

for dy=1:numdyads
    if montageMatch.montagePerm(dy)==1
        for ar=1:numareas
            z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas1(ar)),dy),2); 
            z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas1(ar)),dy),2);
            z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas1(ar)),dy),2); 
            z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas1(ar)),dy),2);

            missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas1(ar)),dy))));
            missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas1(ar)),dy))));
            missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas1(ar)),dy))));
            missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas1(ar)),dy))));
        end 
    elseif montageMatch.montagePerm(dy)==2
        for ar=1:numareas
            z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas2(ar)),dy),2); 
            z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas2(ar)),dy),2);
            z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas2(ar)),dy),2); 
            z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas2(ar)),dy),2);

            missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas2(ar)),dy))));
            missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas2(ar)),dy))));
            missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas2(ar)),dy))));
            missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas2(ar)),dy))));
        end 
    elseif montageMatch.montagePerm(dy)==3
        for ar=1:numareas
            z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas3(ar)),dy),2); 
            z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas3(ar)),dy),2);
            z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas3(ar)),dy),2); 
            z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas3(ar)),dy),2);

            missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas3(ar)),dy))));
            missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas3(ar)),dy))));
            missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas3(ar)),dy))));
            missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas3(ar)),dy))));
        end 
    elseif montageMatch.montagePerm(dy)==4
        if montageMatch.montageOrder(dy)==1
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas1(ar)),dy),2); 
                z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas2(ar)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas1(ar)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas2(ar)),dy),2);

                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas1(ar)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas1(ar)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas2(ar)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas2(ar)),dy))));
            end 
        else
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas2(ar)),dy),2); 
                z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas1(ar)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas2(ar)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas1(ar)),dy),2);

                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas2(ar)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas2(ar)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas1(ar)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas1(ar)),dy))));
            end 
        end
    elseif montageMatch.montagePerm(dy)==5
        if montageMatch.montageOrder(dy)==1
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas2(ar)),dy),2); 
                z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas3(ar)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas2(ar)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas3(ar)),dy),2);

                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas2(ar)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas2(ar)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas3(ar)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas3(ar)),dy))));
            end 
        else
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas3(ar)),dy),2); 
                z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas2(ar)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas3(ar)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas2(ar)),dy),2);

                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas3(ar)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas3(ar)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas2(ar)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas2(ar)),dy))));
            end 
        end
    elseif montageMatch.montagePerm(dy)==6
        if montageMatch.montageOrder(dy)==1
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas1(ar)),dy),2); 
                z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas3(ar)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas1(ar)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas3(ar)),dy),2);

                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas1(ar)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas1(ar)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas3(ar)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas3(ar)),dy))));
            end 
        else
            for ar=1:numareas
                z1_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub1(1:length_diss1,cell2mat(areas3(ar)),dy),2); 
                z2_diss1_areas(:,ar,dy) = nanmean(oxy3D(1).sub2(1:length_diss1,cell2mat(areas1(ar)),dy),2);
                z1_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub1(1:length_diss2,cell2mat(areas3(ar)),dy),2); 
                z2_diss2_areas(:,ar,dy) = nanmean(oxy3D(2).sub2(1:length_diss2,cell2mat(areas1(ar)),dy),2);

                missN_s1(dy,ar,1)=numel(find(isnan(oxy3D(1).sub1(1,cell2mat(areas3(ar)),dy))));
                missN_s1(dy,ar,2)=numel(find(isnan(oxy3D(2).sub1(1,cell2mat(areas3(ar)),dy))));
                missN_s2(dy,ar,1)=numel(find(isnan(oxy3D(1).sub2(1,cell2mat(areas1(ar)),dy))));
                missN_s2(dy,ar,2)=numel(find(isnan(oxy3D(2).sub2(1,cell2mat(areas1(ar)),dy))));
            end 
        end
    end
end