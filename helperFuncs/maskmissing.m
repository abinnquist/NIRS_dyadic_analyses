function [nmask1,nmask2]=maskmissing(montageMatch,missN_s1,missN_s2,numdyads,numareas,areas1,areas2)
loss=0.2;
nmask1=ones(numdyads,numareas,2);
nmask2=ones(numdyads,numareas,2);

for c=1:2
    for dy=1:numdyads
        for a=1:numareas
            if montageMatch.Montage1(dy)==1 
                nChan=length(areas1{a,1});
                cutoff=round(nChan*loss);
                if missN_s1(dy,a,c) > cutoff
                    nmask1(dy,a,c)=0;
                end
            elseif montageMatch.Montage1(dy)==2
                nChan=length(areas2{a,1});
                cutoff=round(nChan*loss);
                if missN_s1(dy,a,c) > cutoff
                    nmask1(dy,a,c)=0;
                end
            end

            if montageMatch.Montage2(dy)==1 
                nChan=length(areas1{a,1});
                cutoff=round(nChan*loss);
                if missN_s2(dy,a,c) > cutoff
                    nmask2(dy,a,c)=0;
                end
            elseif montageMatch.Montage1(dy)==2
                nChan=length(areas2{a,1});
                cutoff=round(nChan*loss);
                if missN_s2(dy,a,c) > cutoff
                    nmask2(dy,a,c)=0;
                end
            end
        end
    end
end

