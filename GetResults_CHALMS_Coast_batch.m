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
MRUNS=30;
%%% ADJUST THIS!!! %%%
ERUNS=1;
batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
    repmat((1:MRUNS)',ERUNS,1)];
batchruns=mat2cell(reshape(1:MRUNS*ERUNS,MRUNS,ERUNS),MRUNS,ones(1,ERUNS));

cd X:\model_results\CHALMS_coast_storm
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('storm',fnamescell(1,:),5);
hind=find(h==1);
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files
load DIST2CBD_east
cd X:\model_results\CHALMS_coast_storm
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
AVGMAP=zeros(NLENGTH,NWIDTH,TMAX,ERUNS);

avgpctdev=zeros(ERUNS,TMAX);
avgthresh=zeros(ERUNS,TMAX);
avgvac=zeros(ERUNS,length(TSTART:TMAX));

vacrlts=zeros(length(hind),length(TSTART:TMAX));
testpctdev=zeros(ERUNS,TMAX);
VARLAYER=zeros(80*80,ERUNS);
for mr=1:length(hind)
    h=strcmp(sprintf('storm_%d_%d.mat',batchind(mr,1),...
        batchind(mr,2)),fnamescell(1,:));
    filename=fnamescell{1,h};
    load(filename)
    
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
    vacrlts(mr,:)=vacrate(TSTART:TMAX);
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
    
%     if strcmp(sprintf('results_CHALMS_Coast_batch_%d_30',batchind(mr,1)),fnamescell(1,hind(mr)))==1
%         buildprob=reshape(VARLAYER,NLENGTH,NWIDTH)/length(hind);
%     end
    VARLAYER(:,batchind(mr,1))=VARLAYER(:,batchind(mr,1))+BASELAYER;    %record frequency of development per cell across model runs
    
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

% Rationale for using Shannon Index to quantify uncertainty in predicting
% housing types for any given location: a given cell could be occupied by 1
% of 8 housing types as well as remain undeveloped. The Shannon index
% describes how rare something is relative to the popularity - should be
% used in combination with a statistic that decribes how consistently a
% cell was occupied in a particualr state.
for ie=1:ERUNS
    ierun=find(batchind(:,1)==ie);
    avgpctdev(ie,:)=mean(pctdev(ierun,:),1);
    avgthresh(ie,:)=avgpctdev(ie,:);
    avgvac(ie,:)=mean(vacrlts(ierun,:),1);
    for tt=TSTART:TMAX
        % probability of housing type across runs, by time step
        htprob_cross(tt,ierun)=mat2cell((cat(2,numltrlts{tt,ierun}).*repmat(z(:,1),1,length(ierun)))./...
            repmat(sum(cat(2,numltrlts{tt,ierun}).*repmat(z(:,1),1,length(ierun)),1),HT,1),HT,ones(1,length(ierun)));
        htentropy_cross(ierun,tt)=-sum(cat(2,htprob_cross{tt,ierun}).*log(cat(2,htprob_cross{tt,ierun})),1)./log(HT);
        
        devreal=cat(2,LTmap{tt,ierun});
        ltprob=histc(devreal,1:HT,2)./length(ierun);
        [maxltprob,imaxtype]=max(ltprob,[],2);
        imaxtype(maxltprob==0)=0;
        
        %%% Fix this to work with LTmap
        devprob=sum((cat(2,LTmap{tt,ierun})~=0),2)./length(ierun);
        
        %%%%% Test representativeness of varmap
        ivarthresh=(devprob >= avgthresh(ie,tt));
        locmat=find(ivarthresh==1);
        testpctdev(ie,tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
        pctdevdiff=testpctdev(ie,tt)-avgpctdev(ie,tt);
        iter=1;
        subpctdevdiff=zeros(1,[]);
        subavgthresh=zeros(1,[]);
        subavgthresh(iter)=avgthresh(ie,tt);
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
            testpctdev(ie,tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
            %         pctdevdiff=testpctdev(tt)-avgpctdev(tt);
            pctdevdiff=testpctdev(ie,tt)-avgpctdev(ie,tt);
            subpctdevdiff(iter+1)=pctdevdiff;
            if abs(subpctdevdiff(iter))-abs(subpctdevdiff(iter+1)) < 0
                if iter > 1
                    ivarthresh=(devprob >= subavgthresh(iter-1));
                    avgthresh(ie,tt)=subavgthresh(iter-1);
                else
                    ivarthresh=(devprob >= avgthresh(ie,tt));
                end
                testpctdev(ie,tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
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
            idevthere=(AVGMAP(:,:,tt-1,ie)~=0);
            indthresh=find(ivarthresh==1);
            inddev=find(idevthere==1);
            idevadd=~ismember(indthresh,inddev);
            subAVGMAP=AVGMAP(:,:,tt-1,ie);
            subAVGMAP(indthresh(idevadd))=imaxtype(indthresh(idevadd));
            AVGMAP(:,:,tt,ie)=subAVGMAP;
        else
            subAVGMAP(ivarthresh)=imaxtype(ivarthresh);
            AVGMAP(:,:,tt,ie)=subAVGMAP;
        end
        testmeandisp=sum(altprop)/length(imapcover);
        %     subCONFMAP=zeros(NLENGTH,NWIDTH);
        %     subCONFMAP(ivarthresh)=confprob(ivarthresh);
        %     CONFMAP(:,:,tt)=subCONFMAP;
    end
end
%% Create figures and results

cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\figs\storm_v1

% resultsfile='baseline_results_070314.txt';
% save(resultsfile,'avgvac','avgpctdev','-ascii')
mdninc=zeros(NLENGTH*NWIDTH,ERUNS);
varinc=zeros(NLENGTH*NWIDTH,ERUNS);
avgbt=zeros(NLENGTH*NWIDTH,ERUNS);
varbt=zeros(NLENGTH*NWIDTH,ERUNS);
avglandsale=zeros(NLENGTH*NWIDTH,ERUNS);
varlandsale=zeros(NLENGTH*NWIDTH,ERUNS);
avgampref=zeros(NLENGTH*NWIDTH,ERUNS);
varampref=zeros(NLENGTH*NWIDTH,ERUNS);
% plot housing types
for j=1:ERUNS
    h1=figure;
    set(h1, 'Color','white','Position',[1,1,700,700],'Visible','off');
    subplot(2,2,1);
    imagesc(AVGMAP(:,:,15,j));
    set(gca,'Visible','off')
    title('t=15');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,2,2);
    imagesc(AVGMAP(:,:,20,j));
    set(gca,'Visible','off')
    title('t=20');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,2,3);
    imagesc(AVGMAP(:,:,25,j));
    set(gca,'Visible','off')
    title('t=25');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,2,4);
    imagesc(AVGMAP(:,:,30,j));
    set(gca,'Visible','off')
    title('t=30');
    set(get(gca,'Title'),'Visible','on')
    saveas(h1,sprintf('avg_housetypes_%d',j),'jpg')
    
    
    % plot consumer incomes
    subdevmap=AVGMAP(:,:,30,j);
    isubdev=find(subdevmap ~= 0);
    for r=1:length(isubdev)
        subincomemap=incomemap(:,batchruns{j});
        inotzero=(subincomemap(isubdev(r),:)~=0);
        mdninc(isubdev(r),j)=median(subincomemap(isubdev(r),inotzero));
        varinc(isubdev(r),j)=var(subincomemap(isubdev(r),inotzero));
    end
    h2=figure;
    orient tall
    set(h2, 'Color','white','OuterPosition',[1,1,400,700],'Visible','off');
    subplot(2,1,1);
    imagesc(reshape(mdninc(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Median Consumer Income');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,1,2);
    imagesc(reshape(varinc(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Variance in Consumer Income');
    set(get(gca,'Title'),'Visible','on')
    saveas(h2,sprintf('final_consumer_income_%d',j),'jpg')
    
    % plot build time
    % isubdev=find(subdevmap ~= 0);
    for r=1:length(isubdev)
        subbtmap=btmap(:,batchruns{j});
        inotzero=(subbtmap(isubdev(r),:)~=0);
        avgbt(isubdev(r),j)=mean(subbtmap(isubdev(r),inotzero));
        varbt(isubdev(r),j)=var(subbtmap(isubdev(r),inotzero));
    end
    h3=figure;
    orient tall
    set(h3, 'Color','white','Position',[1,1,400,700],'Visible','off');
    subplot(2,1,1);
    imagesc(reshape(avgbt(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Average Build Time');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,1,2);
    imagesc(reshape(varbt(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Variance in Build Time');
    set(get(gca,'Title'),'Visible','on')
    saveas(h3,sprintf('build_time_%d',j),'jpg')
    
    % plot entropy of housing type distribution
    h4=figure;
    set(h4, 'Color','white','Position',[1,1,400,400],'Visible','off');
    plot(10:30,mean(htentropy(batchruns{j},10:30),1),'--k','LineWidth',3)
    hold on
    plot(10:30,htentropy(batchruns{j},10:30),'-')
    legend('Mean Entropy','Individual Model Runs')
    xlabel('Time Step')
    ylabel('Entropy of Housing Type Distribution')
    saveas(h4,sprintf('htentropy_%d',j),'jpg')
    % which lot types are leading to each of the trajectories?
    
    % plot land sales
    sublandsales=landsales(:,batchruns{j});
    ilandsale=(sublandsales ~= 0);
    for r=1:NLENGTH*NWIDTH
        avglandsale(r,j)=mean(sublandsales(r,ilandsale(r,:)));
        varlandsale(r,j)=var(sublandsales(r,ilandsale(r,:)));
    end
    h5=figure;
    orient tall
    set(h5, 'Color','white','Position',[700,300,400,700],'Visible','off');
    subplot(2,1,1);
    imagesc(reshape(avglandsale(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Average Land Sales');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,1,2);
    imagesc(reshape(varlandsale(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Variance in Land Sales');
    set(get(gca,'Title'),'Visible','on')
    saveas(h5,sprintf('landsales_%d',j),'jpg')
    
    % plot consumer amenity preferences
    % isubdev=find(subdevmap ~= 0);
    subamprefmap=amprefmap(:,batchruns{j});
    for r=1:length(isubdev)
        inotzero=(subamprefmap(isubdev(r),:)~=0);
        avgampref(isubdev(r),j)=mean(subamprefmap(isubdev(r),inotzero));
        varampref(isubdev(r),j)=var(subamprefmap(isubdev(r),inotzero));
    end
    h6=figure;
    orient tall
    set(h6, 'Color','white','Position',[700,300,400,700],'Visible','off');
    subplot(2,1,1);
    imagesc(reshape(avgampref(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Mean Consumer Amenity Preference');
    set(get(gca,'Title'),'Visible','on')    
    subplot(2,1,2);
    imagesc(reshape(varampref(:,j),NLENGTH,NWIDTH));
    set(gca,'Visible','off')
    colorbar;
    title('Variance in Consumer Amenity Preference');
    set(get(gca,'Title'),'Visible','on')
    saveas(h6,sprintf('amprefmap_%d',j),'jpg')
    
    % plot avg rents
    avgrents=reshape(cell2mat(avgrentstats),HT,TMAX,length(hind));
    mdnrents=median(avgrents(:,:,batchruns{j}),3);
    varrents=var(avgrents(:,:,batchruns{j}),1,3);
    h7=figure;
    set(h7, 'Color','white','Position',[700,300,400,400],'Visible','off');
    plot(10:30,mdnrents(:,10:30),'-');
    title('Median Housing Rents');
    xlabel('Time Step');
    ylabel('Rent');
    legend('HT1','HT2','HT3','HT4','HT5','HT6','HT7','HT8');
    saveas(h7,sprintf('avgrents_%d',j),'jpg')
end
