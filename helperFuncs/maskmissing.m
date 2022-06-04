function [nmask1,nmask2]=maskmissing(montageMatch,missN_s1,missN_s2,numdyads,numareas)
nmask1=ones(numdyads,numareas,2);
nmask2=ones(numdyads,numareas,2);
for convo=1:2
    for dyad=1:numdyads 
        if montageMatch.Montage1(dyad)==1 
            if missN_s1(dyad,1,convo) > 2
                nmask1(dyad,1,convo)=0;
            end
            if missN_s1(dyad,2,convo) > 1
                nmask1(dyad,2,convo)=0;
            end
            if missN_s1(dyad,3,convo) > 3
                nmask1(dyad,3,convo)=0;
            end
            if missN_s1(dyad,4,convo) > 0
                nmask1(dyad,4,convo)=0;
            end
        elseif montageMatch.Montage1(dyad)==2
            if missN_s1(dyad,1,convo) > 1
                nmask1(dyad,1,convo)=0;
            end
            if missN_s1(dyad,2,convo) > 0
                nmask1(dyad,2,convo)=0;
            end
            if missN_s1(dyad,3,convo) > 1
                nmask1(dyad,3,convo)=0;
            end
            if missN_s1(dyad,4,convo) > 1
                nmask1(dyad,4,convo)=0;
            end
        elseif montageMatch.Montage1(dyad)==3
            if missN_s1(dyad,1,convo) > 2
                nmask1(dyad,1,convo)=0;
            end
            if missN_s1(dyad,2,convo) > 1
                nmask1(dyad,2,convo)=0;
            end
            if missN_s1(dyad,3,convo) > 3
                nmask1(dyad,3,convo)=0;
            end
            if missN_s1(dyad,4,convo) > 2
                nmask1(dyad,4,convo)=0;
            end
        end
    end
end

for convo=1:2
    for dyad=1:numdyads 
        if montageMatch.Montage2(dyad)==1 
            if missN_s2(dyad,1,convo) > 2
                nmask2(dyad,1,convo)=0;
            end
            if missN_s2(dyad,2,convo) > 1
                nmask2(dyad,2,convo)=0;
            end
            if missN_s2(dyad,3,convo) > 3
                nmask2(dyad,3,convo)=0;
            end
            if missN_s2(dyad,4,convo) > 0
                nmask2(dyad,4,convo)=0;
            end
        elseif montageMatch.Montage2(dyad)==2
            if missN_s2(dyad,1,convo) > 1
                nmask2(dyad,1,convo)=0;
            end
            if missN_s2(dyad,2,convo) > 0
                nmask2(dyad,2,convo)=0;
            end
            if missN_s2(dyad,3,convo) > 1
                nmask2(dyad,3,convo)=0;
            end
            if missN_s2(dyad,4,convo) > 1
                nmask2(dyad,4,convo)=0;
            end
        elseif montageMatch.Montage2(dyad)==3
            if missN_s2(dyad,1,convo) > 2
                nmask2(dyad,1,convo)=0;
            end
            if missN_s2(dyad,2,convo) > 1
                nmask2(dyad,2,convo)=0;
            end
            if missN_s2(dyad,3,convo) > 3
                nmask2(dyad,3,convo)=0;
            end
            if missN_s2(dyad,4,convo) > 2
                nmask2(dyad,4,convo)=0;
            end
        end
    end
end