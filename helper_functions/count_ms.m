%% summary MS data
function summary_trial=count_ms(MS,tab,path,directionName)
    locationids = 1:8; locationdegrees = {315,135,225,45,270,90,180,0};
    locationlabels = strcat('loc',cellfun(@num2str,locationdegrees,'un',0));
    
    [nTrials, ~] = size(tab);

    a1=[];
    a2=[];
    a3=[];
    a4=[];
    for j = 1 : nTrials
        con=tab(j,9);
        if con == 1 ||  con == 5
            aa=find(MS(:,10)==j);
            a1=[a1;aa];
        elseif con == 2 ||  con == 6
            aa=find(MS(:,10)==j);
            a2=[a2;aa];
        elseif con == 3  ||  con == 7 
            aa=find(MS(:,10)==j);
            a3=[a3;aa];
        else
            aa=find(MS(:,10)==j);
            a4=[a4;aa];
        end
    end
    
    m1=MS(a1,:);
    m2=MS(a2,:);
    m3=MS(a3,:);
    m4=MS(a4,:);
    
    % amplitude
    summary_trial=zeros(10,9);
    summary_trial(1,9)=size(MS,1);
    summary_trial(2,9)=mean(MS(:,11));
    summary_trial(3,9)=max(MS(:,11));
    summary_trial(4,9)=min(MS(:,11));
    % duration
    summary_trial(5,9)=mean(MS(:,2)-MS(:,1));
    summary_trial(6,9)=max(MS(:,2)-MS(:,1));
    summary_trial(7,9)=min(MS(:,2)-MS(:,1));
    % velocity peak
    summary_trial(8,9)=mean(MS(:,3));
    summary_trial(9,9)=max(MS(:,3));
    summary_trial(10,9)=min(MS(:,3));
    
    
    
    con=tab(1,9);
    if con ==1 || con==2 || con==3 || con==4 
        summary_trial(1,1)=size(m1,1);
        %amplitude
        summary_trial(2,1)=mean(m1(:,11));
        summary_trial(3,1)=max(m1(:,11));
        summary_trial(4,1)=min(m1(:,11));
        %duration
        summary_trial(5,1)=mean(m1(:,2)-m1(:,1));
        summary_trial(6,1)=max(m1(:,2)-m1(:,1));
        summary_trial(7,1)=min(m1(:,2)-m1(:,1));
        %velocity
        summary_trial(8,1)=mean(m1(:,3));
        summary_trial(9,1)=max(m1(:,3));
        summary_trial(10,1)=min(m1(:,3));
        
        summary_trial(1,2)=size(m2,1);
        summary_trial(2,2)=mean(m2(:,11));
        summary_trial(3,2)=max(m2(:,11));
        summary_trial(4,2)=min(m2(:,11));
        %duration
        summary_trial(5,2)=mean(m2(:,2)-m2(:,1));
        summary_trial(6,2)=max(m2(:,2)-m2(:,1));
        summary_trial(7,2)=min(m2(:,2)-m2(:,1));
        %velocity
        summary_trial(8,2)=mean(m2(:,3));
        summary_trial(9,2)=max(m2(:,3));
        summary_trial(10,2)=min(m2(:,3));
        
        
        summary_trial(1,3)=size(m3,1);
        summary_trial(2,3)=mean(m3(:,11));
        summary_trial(3,3)=max(m3(:,11));
        summary_trial(4,3)=min(m3(:,11));
        %duration
        summary_trial(5,3)=mean(m3(:,2)-m3(:,1));
        summary_trial(6,3)=max(m3(:,2)-m3(:,1));
        summary_trial(7,3)=min(m3(:,2)-m3(:,1));
        %velocity
        summary_trial(8,3)=mean(m3(:,3));
        summary_trial(9,3)=max(m3(:,3));
        summary_trial(10,3)=min(m3(:,3));
        
        summary_trial(1,4)=size(m4,1);
        summary_trial(2,4)=mean(m4(:,11));
        summary_trial(3,4)=max(m4(:,11));
        summary_trial(4,4)=min(m4(:,11));
         %duration
        summary_trial(5,4)=mean(m4(:,2)-m4(:,1));
        summary_trial(6,4)=max(m4(:,2)-m4(:,1));
        summary_trial(7,4)=min(m4(:,2)-m4(:,1));
        %velocity
        summary_trial(8,4)=mean(m4(:,3));
        summary_trial(9,4)=max(m4(:,3));
        summary_trial(10,4)=min(m4(:,3));
        
    else
        summary_trial(1,5)=size(m1,1);
        summary_trial(2,5)=mean(m1(:,11));
        summary_trial(3,5)=max(m1(:,11));
        summary_trial(4,5)=min(m1(:,11));
        %duration
        summary_trial(5,5)=mean(m1(:,2)-m1(:,1));
        summary_trial(6,5)=max(m1(:,2)-m1(:,1));
        summary_trial(7,5)=min(m1(:,2)-m1(:,1));
        %velocity
        summary_trial(8,5)=mean(m1(:,3));
        summary_trial(9,5)=max(m1(:,3));
        summary_trial(10,5)=min(m1(:,3));
        
        
        summary_trial(1,6)=size(m2,1);
        summary_trial(2,6)=mean(m2(:,11));
        summary_trial(3,6)=max(m2(:,11));
        summary_trial(4,6)=min(m2(:,11));
        %duration
        summary_trial(5,6)=mean(m2(:,2)-m2(:,1));
        summary_trial(6,6)=max(m2(:,2)-m2(:,1));
        summary_trial(7,6)=min(m2(:,2)-m2(:,1));
        %velocity
        summary_trial(8,6)=mean(m2(:,3));
        summary_trial(9,6)=max(m2(:,3));
        summary_trial(10,6)=min(m2(:,3));
        
        
        summary_trial(1,7)=size(m3,1);
        summary_trial(2,7)=mean(m3(:,11));
        summary_trial(3,7)=max(m3(:,11));
        summary_trial(4,7)=min(m3(:,11));
        %duration
        summary_trial(5,7)=mean(m3(:,2)-m3(:,1));
        summary_trial(6,7)=max(m3(:,2)-m3(:,1));
        summary_trial(7,7)=min(m3(:,2)-m3(:,1));
        %velocity
        summary_trial(8,7)=mean(m3(:,3));
        summary_trial(9,7)=max(m3(:,3));
        summary_trial(10,7)=min(m3(:,3));
        
        
        summary_trial(1,8)=size(m4,1);
        summary_trial(2,8)=mean(m4(:,11));
        summary_trial(3,8)=max(m4(:,11));
        summary_trial(4,8)=min(m4(:,11));
        %duration
        summary_trial(5,8)=mean(m4(:,2)-m4(:,1));
        summary_trial(6,8)=max(m4(:,2)-m4(:,1));
        summary_trial(7,8)=min(m4(:,2)-m4(:,1));
        %velocity
        summary_trial(8,8)=mean(m4(:,3));
        summary_trial(9,8)=max(m4(:,3));
        summary_trial(10,8)=min(m4(:,3));
    end
   
end