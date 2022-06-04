function nmask = cfMissingch(missN_s1,missN_s2,numdyads,numareas)
nmask=zeros(numdyads,numareas,2);
for convo=1:2
    if numareas==3
        for dyad=1:numdyads
            for  area=1:numareas
                if area==2 || area==3
                    if missN_s1(dyad,area,convo)>2 || missN_s2(dyad,area,convo)>2
                        nmask(dyad,area,convo)=1;
                    end                
                else
                    if missN_s1(dyad,area,convo)>1 || missN_s2(dyad,area,convo)>1
                        nmask(dyad,area,convo)=1;
                    end
                end
            end
        end
    elseif numareas==5
        for dyad=1:numdyads
            for  area=1:numareas
                if area==1 || area==2
                    if missN_s1(dyad,area,convo)>=1 || missN_s2(dyad,area,convo)>=1
                        nmask(dyad,area,convo)=1;
                    end                
                else
                    if missN_s1(dyad,area,convo)>1 || missN_s2(dyad,area,convo)>1
                        nmask(dyad,area,convo)=1;
                    end
                end
            end
        end
    elseif numareas==6
        for dyad=1:numdyads
            for  area=1:numareas
                if area==1 || area==5
                    if missN_s1(dyad,area,convo)>=1 || missN_s2(dyad,area,convo)>=1
                        nmask(dyad,area,convo)=1;
                    end                
                else
                    if missN_s1(dyad,area,convo)>1 || missN_s2(dyad,area,convo)>1
                        nmask(dyad,area,convo)=1;
                    end
                end
            end
        end
    end
end
