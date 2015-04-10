%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Read-in results files   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This code is written to ingest saved results files from individual
%   model runs, extract targeted statistics, and compile those statistics
%   for visualization and printing to output files (e.g., .csv, .xlsx).

% Reproduce parameters from model so that resutls can be processed
% independently from model runs
NLENGTH=80;
NWIDTH=80;
NCELLS=NLENGTH*NWIDTH;
TSTART=10;
TMAX=30;
HT=8;
cell2mile=0.0395;
cell2ft=cell2mile*5280;
z = 1000*[0.00025    1.5000
    0.00025    2.5000
    0.0005    1.5000
    0.0005    2.5000
    0.0010    1.5000
    0.0010    2.5000
    0.0020    1.5000
    0.0020    2.5000];
% number of model iterations
MRUNS=9;
% number of experimental parameters
NPARMS=1;
% number of experimental variatons on NPARMS
EXPTRUNS=120;
ERUNS=EXPTRUNS;
% index numbers of storm climate settings and 1:MRUNS model runs
batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
    repmat((1:MRUNS)',ERUNS,1)];
batchruns=mat2cell(reshape(1:MRUNS*ERUNS,MRUNS,ERUNS),MRUNS,ones(1,ERUNS));

% Experimental settings
% storm cliamte settings
ma_storm=[0.013 0 0 0 0; 0.0065 0.0065 0 0 0; 0.013 0 0 0 0; 0.0065 0 0 0 0; 0.0584 0.013 0.0065 0 0];
nc_storm=[0.1364 0.0844 0.0714 0.0065 0];
fl_storm=[0.2792 0.2078 0.1753 0.0390 0.0130];
tx_storm=[0.1494 0.1104 0.0779 0.0455 0];

stormprob=cell(4,1);    %[mid-Atlantic average, NC, FL, TX]
stormprob(1)=mat2cell(mean(ma_storm,1),1,5);
stormprob(2)=mat2cell(nc_storm,1,5);
stormprob(3)=mat2cell(fl_storm,1,5);
stormprob(4)=mat2cell(tx_storm,1,5);

batchparms_full=[sum(cat(1,stormprob{1})) sum(cat(1,stormprob{2})) ...
    sum(cat(1,stormprob{3})) sum(cat(1,stormprob{4}))];
batchparms_unq=unique(batchparms_full);
branges=[40000:16000:184000 200001];    % income distribution bins

%%% Adjust this %% 
% navigate to results file storage
cd X:\model_results\CHALMS_coast_gsa_baseline
fnames=dir;
fnamescell=struct2cell(fnames);
% (2) change the prefix of the results file names
h=strncmp('coast_baseline_gsa',fnamescell(1,:),14);
hind=find(h==1);
% add precalculated distance matrix
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files
load DIST2CBD_east
cd X:\model_results\CHALMS_coast_gsa_baseline
distpt=ceil(min(min(dist2cbd)):max(max(dist2cbd)));
density=zeros(NLENGTH*NWIDTH,length(hind));

coastdist=reshape(NWIDTH+1-cumsum(ones(NLENGTH,NWIDTH),2),NCELLS,1);

pctdev=zeros(length(hind),TMAX);
pctdev_zone=zeros(length(hind),TMAX,5);
htprob=cell(TMAX,length(hind));
htentropy=zeros(length(hind),TMAX);
htentropy_zone=zeros(length(hind),TMAX,5);

htprob_cross=cell(TMAX,length(hind));
htentropy_cross=zeros(length(hind),TMAX);

LTmap=cell(TMAX,length(hind));
avgrentstats=cell(1,length(hind));
avgrentstats_zone=zeros(TMAX,length(hind),5);
avgreturnstats=cell(1,length(hind));
avgexptreturnstats=cell(1,length(hind));
numltrlts=cell(TMAX,length(hind));
numltrlts_zone1=cell(TMAX,length(hind));
numltrlts_zone2=cell(TMAX,length(hind));
numltrlts_zone3=cell(TMAX,length(hind));
numltrlts_zone4=cell(TMAX,length(hind));
numltrlts_zone5=cell(TMAX,length(hind));
incomes=cell(1,length(hind));

incomemap=zeros(NLENGTH*NWIDTH,length(hind));
income_zone=zeros(length(hind),5);
amprefmap=zeros(NLENGTH*NWIDTH,length(hind));
housepricemap=zeros(NLENGTH*NWIDTH,length(hind));
btmap=zeros(NLENGTH*NWIDTH,length(hind));
bt_zone=zeros(length(hind),5);
landsales=zeros(NLENGTH*NWIDTH,length(hind));
landsales_zone=zeros(length(hind),TMAX,5);
landsaletime=zeros(NLENGTH*NWIDTH,length(hind));
densitydata=zeros(NLENGTH*NWIDTH,length(hind));
densitydata_zone=zeros(length(hind),TMAX,5);
AVGMAP=zeros(NLENGTH,NWIDTH,TMAX,ERUNS);

avgpctdev=zeros(ERUNS,TMAX);
avgthresh=zeros(ERUNS,TMAX);
avgvac=zeros(ERUNS,length(TSTART:TMAX));

rentrlts=cell(length(hind));
totdam=zeros(length(hind),TMAX);
totdam_zone=zeros(length(hind),TMAX,5);

vacrlts=zeros(length(hind),length(TSTART:TMAX));
vacrlts_zone=zeros(length(hind),5);
testpctdev=zeros(ERUNS,TMAX);
VARLAYER=zeros(80*80,ERUNS);
indmap=reshape(1:NLENGTH*NWIDTH,NLENGTH,NWIDTH);
izone1=indmap(:,1:50);
izone3=indmap(31:NLENGTH,51:70);
for mr=1:length(hind)   % MRUNS*EXPTRUNS
    h=strcmp(sprintf('coast_baseline_gsa%d_%d.mat',batchind(mr,1),...
        batchind(mr,2)),fnamescell(1,:));
    filename=fnamescell{1,h};
    load(filename)
    
    damagemap=0.01.*(X{mr,2}*10.23749-0.23462*(coastdist*cell2ft/1000)+...
            0.001649*(coastdist*cell2ft/1000).^2);
    
    izone5=find(BUILDTIME == TSTART);
    
    zone2inds=indmap(1:30,51:70);
    izone2=zone2inds(~ismember(zone2inds,izone5));
    zone4inds=indmap(:,71:80);
    izone4=zone4inds(~ismember(zone4inds,izone5));
    
    % lot type map over time
    LTmap(:,mr)=mat2cell(LOTTYPE,NLENGTH*NWIDTH,ones(1,TMAX));
    devarea=(cat(2,LTmap{:,mr})~=0);
    pctdev(mr,:)=sum(devarea,1)/(NLENGTH*NWIDTH);
    
    % probability of housing type being built in a single run
    numltrlts(:,mr)=mat2cell(numlt,HT,ones(1,TMAX));
    htprob(:,mr)=mat2cell((numlt.*repmat(z(:,1),1,TMAX))./...
        repmat(sum(numlt.*repmat(z(:,1),1,TMAX),1),HT,1),HT,ones(1,TMAX));
    htentropy(mr,:)=-sum(cat(2,htprob{:,mr}).*log(cat(2,htprob{:,mr})),1)./log(HT);
    
    % find density gradient
    subltmap=cat(2,LTmap{:,mr});
    densitydata(:,mr)=subltmap(:,TMAX);
    
    for d=1:length(distpt)-1
        lotsize=densitydata(find(dist2cbd >= distpt(d) & dist2cbd < distpt(d+1)),mr);
        inotzero=(lotsize~=0);
        if isempty(find(inotzero,1))==1
            density(d,mr)=0;
        else
            density(d,mr)=1/(sum(1./z(lotsize(inotzero),1))/length(lotsize(inotzero)));
        end
    end

    % time series data
    vacrlts(mr,:)=vacrate(TSTART:TMAX); % vacancy rate
    avgrentstats(mr)=mat2cell(avgrent,HT,TMAX); % average rents
    Realavgret(isnan(Realavgret))=0;   % average return per housetype based on sale price
    avgreturnstats(mr)=mat2cell(Realavgret,HT,TMAX);
    Realavgexptret(isnan(Realavgexptret))=0;    % average expected return based on asking price
    avgexptreturnstats(mr)=mat2cell(Realavgexptret,HT,TMAX);
    rentrlts(mr)=mat2cell(RENT(:,TSTART:TMAX),NCELLS,length(TSTART:TMAX));
    % consumer incomes at TMAX
    incomes(mr)=mat2cell(cat(1,CONINFO{:,1}),length(CONINFO),1);
    
    % mapped variables
    btmap(:,mr)=BUILDTIME;
    inotvac=find(cat(1,lotchoice{:,4})==1);
    for i=1:length(cat(1,Lottype{inotvac,1}))
        incomemap(Lottype{inotvac(i),2},mr)=CONINFO{lotchoice{inotvac(i),5},1}*...
            ones(length(Lottype{inotvac(i),2}),1);  % incomes
        amprefmap(Lottype{inotvac(i),2},mr)=CONINFO{lotchoice{inotvac(i),5},6}*...
            ones(length(Lottype{inotvac(i),2}),1);  % consumer preference for amenity
        housepricemap(Lottype{inotvac(i),2},mr)=lotchoice{inotvac(i),7};    %house asking/sale prices
    end
    
    % land sales
    for nf=1:length(farmsoldinfo(:,1))
        ifarm=find(cat(1,LANDINFO{1,10})==farmsoldinfo(nf,1));
        landsales(ifarm,mr)=farmsoldinfo(nf,4);
        landsaletime(ifarm,mr)=farmsoldinfo(nf,2);
    end
    
    for tt=TSTART:TMAX
        totdam(mr,tt)=sum(RENT(:,tt).*damagemap);
    end
    
    for iz=1:5
        zonemap=zeros(NLENGTH,NWIDTH);
        for it=1:TMAX
            if iz==1
                zonemap(izone1)=LOTTYPE(izone1,it);
                idevarea=(zonemap~=0);
                islot=ismember(lotchoice(:,2),find(idevarea==1));
                numltrlts_zone1(it,mr)=mat2cell(histc(zonemap(idevarea),1:HT),HT,1);
                htentropy_zone(mr,it,iz)=-sum(((numltrlts_zone1{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone1{it,mr}.*z(:,1))).*log(((numltrlts_zone1{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone1{it,mr}.*z(:,1)))),1)./log(HT);
                zsize=size(izone1);
                pctdev_zone(mr,it,iz)=length(find(idevarea==1))/(zsize(1)*zsize(2));
            elseif iz==2
                zonemap(izone2)=LOTTYPE(izone2,it);
                idevarea=(zonemap~=0);
                islot=ismember(lotchoice(:,2),find(idevarea==1));
                numltrlts_zone2(it,mr)=mat2cell(histc(zonemap(idevarea),1:HT),HT,1);
                htentropy_zone(mr,it,iz)=-sum(((numltrlts_zone2{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone2{it,mr}.*z(:,1))).*log(((numltrlts_zone2{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone2{it,mr}.*z(:,1)))),1)./log(HT);
                zsize=size(izone2);
                pctdev_zone(mr,it,iz)=length(find(idevarea==1))/(zsize(1)*zsize(2));
            elseif iz==3
                zonemap(izone3)=LOTTYPE(izone3,it);
                idevarea=(zonemap~=0);
                islot=ismember(lotchoice(:,2),find(idevarea==1));
                numltrlts_zone3(it,mr)=mat2cell(histc(zonemap(idevarea),1:HT),HT,1);
                htentropy_zone(mr,it,iz)=-sum(((numltrlts_zone3{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone3{it,mr}.*z(:,1))).*log(((numltrlts_zone3{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone3{it,mr}.*z(:,1)))),1)./log(HT);
                zsize=size(izone3);
                pctdev_zone(mr,it,iz)=length(find(idevarea==1))/(zsize(1)*zsize(2));
            elseif iz==4
                zonemap(izone4)=LOTTYPE(izone4,it);
                idevarea=(zonemap~=0);
                islot=ismember(lotchoice(:,2),find(idevarea==1));
                numltrlts_zone4(it,mr)=mat2cell(histc(zonemap(idevarea),1:HT),HT,1);
                htentropy_zone(mr,it,iz)=-sum(((numltrlts_zone4{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone4{it,mr}.*z(:,1))).*log(((numltrlts_zone4{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone4{it,mr}.*z(:,1)))),1)./log(HT);
                zsize=size(izone4);
                pctdev_zone(mr,it,iz)=length(find(idevarea==1))/(zsize(1)*zsize(2));
            elseif iz==5
                zonemap(izone5)=LOTTYPE(izone5,it);
                idevarea=(zonemap~=0);
                islot=ismember(lotchoice(:,2),find(idevarea==1));
                numltrlts_zone5(it,mr)=mat2cell(histc(zonemap(idevarea),1:HT),HT,1);
                htentropy_zone(mr,it,iz)=-sum(((numltrlts_zone5{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone5{it,mr}.*z(:,1))).*log(((numltrlts_zone5{it,mr}.*z(:,1))./...
                    sum(numltrlts_zone5{it,mr}.*z(:,1)))),1)./log(HT);
                zsize=size(izone5);
                pctdev_zone(mr,it,iz)=length(find(idevarea==1))/(zsize(1)*zsize(2));
            end
            densitydata_zone(mr,it,iz)=mean(z(LOTTYPE(idevarea,it),1));
            idev=find(idevarea==1);
            ilsale=ismember(idevarea,find(landsaletime(idevarea)==it));
            landsales_zone(mr,it,iz)=mean(landsales(idevarea(ilsale),mr));
            totdam_zone(mr,it,iz)=sum(RENT(idevarea,it).*damagemap(idevarea,it));
        end
        income_zone(mr,iz)=median(incomemap(idevarea,mr));
        avgrentstats_zone(:,mr,iz)=mean(lotchoice(ismember(lotchoice(:,2),idevarea),7));
        bt_zone(mr,iz)=mean(BUILDTIME(idevarea));
    end
    
    % build development proability map
    VARLAYER(:,batchind(mr,1))=VARLAYER(:,batchind(mr,1))+BASELAYER;    %record frequency of development per cell across model runs

    % clear single model run results
    clear('consumerstats','vacstats','BUILDTIME','VACLAND','RENT','RETURN',...
        'LOTTYPE','BASELAYER','Rpop','Rvacrate','Rvaclots',...
        'numlt','Rleftoverpop','avgrentdynms','rentdynstats',...
        'farmsoldinfo','avgrent','avgfarmsize','stdfarmsize','DELTA',...
        'survivalrate','LOCWGHT','REGWGHT','PCTSEARCH','zeta','HIBETA',...
        'MIDBETA','LOWBETA','POPGROW','ccost','newhouseset','bidtot',...
        'meandisp','maxdevdist','setupmap','zonedensity',...
        'vacrate','Farminfo','oldincome','Realreturn',...
        'Realavgret','Exptrentdiff','Avgexptdiff','htincome',...
        'numtotbids','totfarmrecord','htperyear','Newbidlevel',...
        'totbrokerrecord','Farmdist2dev','Bprojmap','Bmodelmap',...
        'WTAlandmap','WTPlandmap','Landprojmap','Landmodmap',...
        'bidshare','buildshare','landdemand','EXPTHOUSE','lotchoice',...
        'Exptprofit','Exptret','Realexptret','Realavgexptret','idealset',...
        'profset','avgbrokervar','carrycost','Lottype','CONINFO','PREFMAP',...
        'TSI','IMPACT','DAMAGE','LANDINFO','lotlocate')
end
save baseline_gsa_results avgrentstats avgrentstats_zone bt_zone btmap ...
    densitydata densitydata_zone htentropy htentropy_zone income_zone ...
    incomemap landsales landsales_zone numltrlts numltrlts_zone1 ...
    numltrlts_zone2 numltrlts_zone3 numltrlts_zone4 numltrlts_zone5 ...
    pctdev pctdev_zone totdam totdam_zone vacrlts

