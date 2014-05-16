%%%%%%% Data input %%%%%%%%%%%%%%%

resnum=[3963; 10230; 7892; 3947; 2096; 1255];
rtstart(:,1)=[1; 3; 8; 18; 28; 38; 48];
rtstart(:,2)=[2; 7; 17; 27; 37; 47; 50];
restimedata=zeros([],1);
for i=1:length(resnum)
    restimedata(length(restimedata)+1:length(restimedata)+resnum(i),1)=...
        rtstart(i,1)+(rtstart(i,2)-rtstart(i,1))*rand(resnum(i),1);
end
avgrestime=mean(restimedata);
stdrestime=std(restimedata);

%%% construction costs %%%
% costset=[201624.1197
% 264124.1197
% 326624.1197
% 208273.7676
% 270773.7676
% 333273.7676
% 217748.2394
% 280248.2394
% 342748.2394
% 232872.3592
% 295372.3592
% 357872.3592
% 238120.5986
% 300620.5986
% 363120.5986
% 256368.838
% 318868.838
% 381368.838];
inflate=207.342/113.6; %inflation adjustment from 1987 to 2007 dollars
% presqftcost=linspace(85,165,HT); %mean cost of 125
presqftcost=100*ones(1,HT); %mean cost of 125
streetcost=[7000; 7000; 7000; 15000; 15000; 15000; 25000; 25000; 25000]; %1987 dollars
sewercost=[8000; 8000; 8000; 18000; 18000; 18000; 5000; 5000; 5000];   %2007 dollars
stepsize=length(presqftcost)/(HT);
stepcount=0;
sqftcost=zeros(HT,1);
for hh=1:HT
    sqftcost(hh,1)=mean(presqftcost(stepcount+1:stepcount+stepsize));
    strtcst(hh,1)=mean(streetcost(stepcount+1:stepcount+stepsize));
    swrcst(hh,1)=mean(sewercost(stepcount+1:stepcount+stepsize));
    stepcount=stepcount+stepsize;
end

infracost(1:HT,1)=strtcst.*inflate+swrcst;

%between 10000:49000, 50000:99999, 100000:greater
incomenum=round(116783.*[0.50 .35 0.15]);
inspan(:,1)=[40000; 80000; 120000];
inspan(:,2)=[69999; 119999; 200000];
% incomenum=round(116783.*[0.405 .393 0.202]);
% inspan(:,1)=[20000; 39334; 100615];
% inspan(:,2)=[39333; 100614; 200000];
incomedata=zeros([],1);
for i=1:length(incomenum)
    incomedata(length(incomedata)+1:length(incomedata)+incomenum(i),1)=...
        inspan(i,1)+(inspan(i,2)-inspan(i,1))*rand(incomenum(i),1);
end
% muhat=expfit(incomedata);
parmhat=lognfit(log(incomedata));
avgincome=mean(incomedata);
stdincome=std(incomedata);
% avgincome=101376;
% stdincome=3780;

% betanum=[7240; 3959; 2939; 2069; 4358];
% bspan(:,1)=[0.1; 0.2; 0.25; 0.3; 0.35];
% bspan(:,2)=[0.19; 0.249; 0.299; 0.349; 0.45];
% betadata=zeros([],1);
% for i=1:length(betanum)
%     betadata(length(betadata)+1:length(betadata)+betanum(i),1)=...
%         bspan(i,1)+(bspan(i,2)-bspan(i,1))*rand(betanum(i),1);
% end

%%%%%%%%%% Farmers %%%%%%%%%%%%
% From USDA Census data for MD, PA, VA, DE
returnmeans=[201.17 150.17 51.82 570.94];
sizemans=[52 65 70 31];
AVGFARMRETURN=2486.3;
STDFARMRETURN=AVGFARMRETURN*0.10;