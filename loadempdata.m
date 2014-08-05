function [restimedata,avgrestime,stdrestime,strtcst,swrcst,infracost,...
    incomedata,parmhat,avgincome,stdincome]=loadempdata(HT,rtstart,resnum,...
    presqftcost,streetcost,sewercost,inflate,incomenum,inspan)

restimedata=zeros([],1);
for i=1:length(resnum)
    restimedata(length(restimedata)+1:length(restimedata)+resnum(i),1)=...
        rtstart(i,1)+(rtstart(i,2)-rtstart(i,1))*rand(resnum(i),1);
end
avgrestime=mean(restimedata);
stdrestime=std(restimedata);

stepsize=length(presqftcost)/(HT);
stepcount=0;
sqftcost=zeros(HT,1);
strtcst=zeros(HT,1);
swrcst=zeros(HT,1);
for hh=1:HT
    sqftcost(hh,1)=mean(presqftcost(stepcount+1:stepcount+stepsize));
    strtcst(hh,1)=mean(streetcost(stepcount+1:stepcount+stepsize));
    swrcst(hh,1)=mean(sewercost(stepcount+1:stepcount+stepsize));
    stepcount=stepcount+stepsize;
end
infracost(1:HT,1)=strtcst.*inflate+swrcst;

incomedata=zeros([],1);
for i=1:length(incomenum)
    incomedata(length(incomedata)+1:length(incomedata)+incomenum(i),1)=...
        inspan(i,1)+(inspan(i,2)-inspan(i,1))*rand(incomenum(i),1);
end
parmhat=lognfit(log(incomedata));
avgincome=mean(incomedata);
stdincome=std(incomedata);

end