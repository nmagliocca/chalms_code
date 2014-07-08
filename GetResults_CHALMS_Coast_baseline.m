%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Read-in results files   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NLENGTH=80;
NWIDTH=80;
TSTART=10;
TMAX=30;
HT=8;
z = 1000*[0.00025    1.5000
    0.00025    2.5000
    0.0005    1.5000
    0.0005    2.5000
    0.0010    1.5000
    0.0010    2.5000
    0.0020    1.5000
    0.0020    2.5000];

cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\results
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('results_CHALMS_Coast_baseline',fnamescell(1,:),29);
hind=find(h==1);
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files
load DIST2CBD_east
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\results
distpt=ceil(min(min(dist2cbd)):max(max(dist2cbd)));
density=zeros(NLENGTH*NWIDTH,length(hind));

% buildprob=zeros(NLENGTH,NWIDTH);
pctdev=zeros(length(hind),TMAX);
htprob=cell(TMAX,length(hind));
htentropy=zeros(length(hind),TMAX);

htprob_cross=cell(TMAX,length(hind));
htentropy_cross=zeros(length(hind),TMAX);

LTmap=cell(TMAX,length(hind));
avgrentstats=cell(1,length(hind));
numltrlts=cell(TMAX,length(hind));

incomemap=zeros(NLENGTH*NWIDTH,length(hind));
amprefmap=zeros(NLENGTH*NWIDTH,length(hind));
housepricemap=zeros(NLENGTH*NWIDTH,length(hind));
btmap=zeros(NLENGTH*NWIDTH,length(hind));
landsales=zeros(NLENGTH*NWIDTH,length(hind));
densitydata=zeros(NLENGTH*NWIDTH,length(hind));
AVGMAP=zeros(NLENGTH,NWIDTH,TMAX);

vacrlts=zeros(length(TSTART:TMAX),length(hind));
testpctdev=zeros(1,TMAX);
for mr=1:length(hind)
    filename=char(fnamescell(1,hind(mr)));
    load(eval('filename'))
    
    LTmap(:,mr)=mat2cell(LOTTYPE,NLENGTH*NWIDTH,ones(1,TMAX));
    devarea=(cat(2,LTmap{:,mr})~=0);
    pctdev(mr,:)=sum(devarea,1)/(NLENGTH*NWIDTH);
    
    % probability of housing type in within a single run
    numltrlts(:,mr)=mat2cell(numlt,HT,ones(1,TMAX));
    htprob(:,mr)=mat2cell((numlt.*repmat(z(:,1),1,TMAX))./...
        repmat(sum(numlt.*repmat(z(:,1),1,TMAX),1),HT,1),HT,ones(1,TMAX));
%     fullset=[1-pctdev(mr,:); repmat(pctdev(mr,:),HT,1).*cat(2,htprob{:,mr})];
%     entropy(:,mr)=-sum(fullset.*log(fullset),1)./...
%         log(HT+1);
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
    vacrlts(:,mr)=vacrate(TSTART:TMAX);
    avgrentstats(mr)=mat2cell(avgrent,HT,TMAX);
    
    % mapped variables
    btmap(:,mr)=BUILDTIME;
    inotvac=find(cat(1,lotchoice{:,4})==1);
    for i=1:length(cat(1,Lottype{inotvac,1}))
        incomemap(Lottype{inotvac(i),2},mr)=CONINFO{lotchoice{inotvac(i),5},1}*...
            ones(length(Lottype{inotvac(i),2}),1);
        amprefmap(Lottype{inotvac(i),2},mr)=CONINFO{lotchoice{inotvac(i),5},6}*...
            ones(length(Lottype{inotvac(i),2}),1);
        housepricemap(Lottype{inotvac(i),2},mr)=lotchoice{inotvac(i),7};
    end
    
    
    for nf=1:length(farmsoldinfo(:,1))
        ifarm=find(cat(1,LANDINFO{1,10})==farmsoldinfo(nf,1));
        landsales(ifarm,mr)=farmsoldinfo(nf,4);
    end
    
    if strncmp('results_CHALMS_Coast_baseline_30',fnamescell(1,hind(mr)),32)==1
        buildprob=reshape(VARLAYER,NLENGTH,NWIDTH)/length(hind);
    end

    clear('consumerstats','vacstats','BUILDTIME','VACLAND','RENT','RETURN',...
        'LOTTYPE','BASELAYER','Rpop','Rvacrate','Rvaclots',...
        'numlt','Rleftoverpop','avgrentdynms','rentdynstats',...
        'farmsoldinfo','avgrent','avgfarmsize','stdfarmsize','DELTA',...
        'survivalrate','LOCWGHT','REGWGHT','PCTSEARCH','zeta','HIBETA',...
        'MIDBETA','LOWBETA','POPGROW','ccost','newhouseset','bidtot',...
        'meandisp','maxdevdist','setupmap','zonedensity',...
        'VARLAYER','vacrate','Farminfo','oldincome','Realreturn',...
        'Realavgret','Exptrentdiff','Avgexptdiff','htincome',...
        'numtotbids','totfarmrecord','htperyear','Newbidlevel',...
        'totbrokerrecord','Farmdist2dev','Bprojmap','Bmodelmap',...
        'WTAlandmap','WTPlandmap','Landprojmap','Landmodmap',...
        'bidshare','buildshare','landdemand','EXPTHOUSE','lotchoice',...
        'Exptprofit','Exptret','Realexptret','Realavgexptret','idealset',...
        'profset','avgbrokervar','carrycost','Lottype','CONINFO','PREFMAP',...
        'TSI','IMPACT','DAMAGE','LANDINFO','lotlocate')
end
% buildprob=reshape(VARLAYER,NLENGTH,NWIDTH)/length(hind);

avgpctdev=mean(pctdev,1);
avgthresh=avgpctdev;
% Rationale for using Shannon Index to quantify uncertainty in predicting
% housing types for any given location: a given cell could be occupied by 1
% of 8 housing types as well as remain undeveloped. The Shannon index
% describes how rare something is relative to the popularity - should be
% used in combination with a statistic that decribes how consistently a
% cell was occupied in a particualr state.
for tt=TSTART:TMAX
    % probability of housing type across runs, by time step
    htprob_cross(tt,:)=mat2cell((cat(2,numltrlts{tt,:}).*repmat(z(:,1),1,length(hind)))./...
        repmat(sum(cat(2,numltrlts{tt,:}).*repmat(z(:,1),1,length(hind)),1),HT,1),HT,ones(1,length(hind)));
    htentropy_cross(mr,:)=-sum(cat(2,htprob_cross{tt,:}).*log(cat(2,htprob_cross{tt,:})),1)./log(HT);
    
    devreal=cat(2,LTmap{tt,:});
    ltprob=histc(devreal,1:HT,2)./length(hind);
    [maxltprob,imaxtype]=max(ltprob,[],2);
    imaxtype(maxltprob==0)=0;
    
    %%% Fix this to work with LTmap
    devprob=sum((cat(2,LTmap{tt,:})~=0),2)./length(hind);

    %%%%% Test representativeness of varmap
    ivarthresh=(devprob >= avgthresh(tt));
    locmat=find(ivarthresh==1);
    testpctdev(tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
    pctdevdiff=testpctdev(tt)-avgpctdev(tt);
    iter=1;
    subpctdevdiff=zeros(1,[]);
    subavgthresh=zeros(1,[]);
    subavgthresh(iter)=avgthresh(tt);
    while abs(pctdevdiff) > 0.01
        subpctdevdiff(iter)=pctdevdiff;
        if iter > 1
            subavgthresh(iter)=subavgthresh(iter-1)+subpctdevdiff(iter);
        else
            subavgthresh(iter)=avgthresh(tt)+subpctdevdiff(iter);
        end
%         avgthresh(tt)=avgthresh(tt)+pctdevdiff;
        
%         ivarthresh=(devprob >= avgthresh(tt));
        ivarthresh=(devprob >= subavgthresh(iter));
        locmat=find(ivarthresh==1);
        testpctdev(tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
%         pctdevdiff=testpctdev(tt)-avgpctdev(tt);
        pctdevdiff=testpctdev(tt)-avgpctdev(tt);
        subpctdevdiff(iter+1)=pctdevdiff;
        if abs(subpctdevdiff(iter))-abs(subpctdevdiff(iter+1)) < 0
            if iter > 1
                ivarthresh=(devprob >= subavgthresh(iter-1));
                avgthresh(tt)=subavgthresh(iter-1);
            else
                ivarthresh=(devprob >= avgthresh(tt));
            end
            testpctdev(tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
            break
        end
        iter=iter+1;
        %%% Need to include a kill switch if abs(pctdevdiff) is not getting
        %%% smaller - need temporal tracking of pctdevdiff with (iter)
        %%% Need to track avgthresh so that previous state can be restored
        %%% if performance does not improve
    end

    imapcover=find(ivarthresh==1);
    testmap=zeros(NLENGTH,NWIDTH);
    testmap(imapcover)=1;
    altprop=zeros(length(find(imapcover)),1);
    for ic=1:length(imapcover)
        [row,col]=ind2sub([NLENGTH NWIDTH],imapcover(ic));
        updir=max(row-2,1);
        dndir=min(row+2,NLENGTH);
        lfdir=max(col-2,1);
        rtdir=min(col+2,NWIDTH);
        altprop(ic)=length(find(testmap(updir:dndir,lfdir:rtdir)==0))/...
            (length(updir:dndir)*length(lfdir:rtdir));
    end
    subAVGMAP=zeros(NLENGTH,NWIDTH);
    if tt > TSTART
        idevthere=(AVGMAP(:,:,tt-1)~=0);
        indthresh=find(ivarthresh==1);
        inddev=find(idevthere==1);
        idevadd=~ismember(indthresh,inddev);
        subAVGMAP=AVGMAP(:,:,tt-1);
        subAVGMAP(indthresh(idevadd))=imaxtype(indthresh(idevadd));
        AVGMAP(:,:,tt)=subAVGMAP;
    else
        subAVGMAP(ivarthresh)=imaxtype(ivarthresh);
        AVGMAP(:,:,tt)=subAVGMAP;
    end
    testmeandisp=sum(altprop)/length(imapcover);
%     subCONFMAP=zeros(NLENGTH,NWIDTH);
%     subCONFMAP(ivarthresh)=confprob(ivarthresh);
%     CONFMAP(:,:,tt)=subCONFMAP;
end

%% Create figures and results

avgvac=mean(vacrlts,1);

resultsfile='baseline_results_070314.txt';
save(resultsfile,'avgvac','avgpctdev','-ascii')

% plot housing types
figure(1)
set(gcf, 'Color','white')
subplot(2,2,1)
imagesc(AVGMAP(:,:,15))
title('t=15')
subplot(2,2,2)
imagesc(AVGMAP(:,:,20))
title('t=20')
subplot(2,2,3)
imagesc(AVGMAP(:,:,25))
title('t=25')
subplot(2,2,4)
imagesc(AVGMAP(:,:,30))
title('t=30')

% plot consumer incomes
subdevmap=AVGMAP(:,:,30);
mdninc=zeros(NLENGTH*NWIDTH,1);
varinc=zeros(NLENGTH*NWIDTH,1);
isubdev=find(subdevmap ~= 0);
for r=1:length(isubdev)
    inotzero=(incomemap(isubdev(r),:)~=0);
    mdninc(isubdev(r))=median(incomemap(isubdev(r),inotzero));
    varinc(isubdev(r))=var(incomemap(isubdev(r),inotzero));
end
figure(2)
set(gcf, 'Color','white')
subplot(2,1,1)
imagesc(reshape(mdninc,NLENGTH,NWIDTH))
colorbar
title('Median Consumer Income')
subplot(2,1,2)
imagesc(reshape(varinc,NLENGTH,NWIDTH))
colorbar
title('Variance in Consumer Income')

% plot build time
avgbt=zeros(NLENGTH*NWIDTH,1);
varbt=zeros(NLENGTH*NWIDTH,1);
isubdev=find(subdevmap ~= 0);
for r=1:length(isubdev)
    inotzero=(btmap(isubdev(r),:)~=0);
    avgbt(isubdev(r))=mean(btmap(isubdev(r),inotzero));
    varbt(isubdev(r))=var(btmap(isubdev(r),inotzero));
end
figure(3)
set(gcf, 'Color','white')
subplot(2,1,1)
imagesc(reshape(avgbt,NLENGTH,NWIDTH))
colorbar
title('Average Build Time')
subplot(2,1,2)
imagesc(reshape(varbt,NLENGTH,NWIDTH))
colorbar
title('Variance in Build Time')

% plot entropy of housing type distribution
figure(4)
set(gcf, 'Color','white')
plot(10:30,mean(htentropy(:,10:30),1),'--k','LineWidth',3)
hold on
plot(10:30,htentropy(:,10:30),'-')
legend('Mean Entropy','Individual Model Runs')
xlabel('Time Step')
ylabel('Entropy of Housing Type Distribution')
% which lot types are leading to each of the trajectories?

% plot land sales
avglandsale=zeros(NLENGTH*NWIDTH,1);
varlandsale=zeros(NLENGTH*NWIDTH,1);
ilandsale=(landsales ~= 0);
for r=1:NLENGTH*NWIDTH
    avglandsale(r)=mean(landsales(r,ilandsale(r,:)));
    varlandsale(r)=var(landsales(r,ilandsale(r,:)));
end
figure(5)
set(gcf, 'Color','white')
subplot(2,1,1)
imagesc(reshape(avglandsale,NLENGTH,NWIDTH))
colorbar
title('Average Land Sales')
subplot(2,1,2)
imagesc(reshape(varlandsale,NLENGTH,NWIDTH))
colorbar
title('Variance in Land Sales')

% plot consumer amenity preferences
avgampref=zeros(NLENGTH*NWIDTH,1);
varampref=zeros(NLENGTH*NWIDTH,1);
isubdev=find(subdevmap ~= 0);
for r=1:length(isubdev)
    inotzero=(amprefmap(isubdev(r),:)~=0);
    avgampref(isubdev(r))=mean(amprefmap(isubdev(r),inotzero));
    varampref(isubdev(r))=var(amprefmap(isubdev(r),inotzero));
end
figure(6)
set(gcf, 'Color','white')
subplot(2,1,1)
imagesc(reshape(avgampref,NLENGTH,NWIDTH))
colorbar
title('Mean Consumer Amenity Preference')
subplot(2,1,2)
imagesc(reshape(varampref,NLENGTH,NWIDTH))
colorbar
title('Variance in Consumer Amenity Preference')

% plot avg rents
avgrents=reshape(cell2mat(avgrentstats),HT,TMAX,length(hind));
mdnrents=median(avgrents,3);
varrents=var(avgrents,1,3);
figure(7)
set(gcf, 'Color','white')
plot(10:30,mdnrents(:,10:30),'-')
title('Median Housing Rents')
xlabel('Time Step')
ylabel('Rent')
    
