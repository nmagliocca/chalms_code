%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Read-in results files   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NLENGTH=80;
NWIDTH=80;
NCELLS=NLENGTH*NWIDTH;
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
NPARMS=1;
EXPTRUNS=4;
ERUNS=EXPTRUNS^NPARMS;
batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
    repmat((1:MRUNS)',ERUNS,1)];
batchruns=mat2cell(reshape(1:MRUNS*ERUNS,MRUNS,ERUNS),MRUNS,ones(1,ERUNS));

% Experimental Parms
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
branges=[40000:16000:184000 200001];
% baseruns=511:540;

% matguide=batchparms_full;
% for ii=1:EXPTRUNS
% matguide(batchparms_full(:,1)==batchparms_unq(ii,1),1)=ii;
% matguide(batchparms_full(:,2)==batchparms_unq(ii,2),2)=ii;
% end

cd X:\model_results\CHALMS_alt_storm_climate
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('storm_clim',fnamescell(1,:),10);
hind=find(h==1);
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files
load DIST2CBD_east
cd X:\model_results\CHALMS_alt_storm_climate
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
incomes=cell(1,length(hind));

incomemap=zeros(NLENGTH*NWIDTH,length(hind));
amprefmap=zeros(NLENGTH*NWIDTH,length(hind));
housepricemap=zeros(NLENGTH*NWIDTH,length(hind));
btmap=zeros(NLENGTH*NWIDTH,length(hind));
landsales=zeros(NLENGTH*NWIDTH,length(hind));
landsaletime=zeros(NLENGTH*NWIDTH,length(hind));
densitydata=zeros(NLENGTH*NWIDTH,length(hind));
AVGMAP=zeros(NLENGTH,NWIDTH,TMAX,ERUNS);

avgpctdev=zeros(ERUNS,TMAX);
avgthresh=zeros(ERUNS,TMAX);
avgvac=zeros(ERUNS,length(TSTART:TMAX));

vacrlts=zeros(length(hind),length(TSTART:TMAX));
testpctdev=zeros(ERUNS,TMAX);
VARLAYER=zeros(80*80,ERUNS);
for mr=1:length(hind)
    h=strcmp(sprintf('storm_clim%d_%d.mat',batchind(mr,1),...
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
    
    incomes(mr)=mat2cell(cat(1,CONINFO{:,1}),length(CONINFO),1);
    
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
        landsaletime(ifarm,mr)=farmsoldinfo(nf,2);
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
                subavgthresh(iter)=avgthresh(ie,tt)+subpctdevdiff(iter);
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

%% Single run comparisons
irun1=(batchind(:,2)==1);
ls1=landsales(:,batchind(:,2)==1);
bt1=btmap(:,batchind(:,2)==1);
inc1=incomemap(:,batchind(:,2)==1);

inc_t=zeros(10,TMAX,4);
avgrent_t=zeros(HT,TMAX,4);
vacrates=zeros(4,21);
landprice=zeros(NCELLS,4);
landtime=zeros(NCELLS,4);
landsale_run=cell(1,4);
saletime=cell(1,4);
branges=[40000:16000:184000 200001];
htype_run=zeros(HT,TMAX,4);
rentmat=zeros(4,TMAX,HT);
inc_wghavg=zeros(TMAX,4);
for i=1:4
    for it=TSTART:TMAX
        ibt=(bt1(:,i)==it);
        occinc=unique(inc1(ibt,i));
        occinc=occinc(occinc~=0);
        h=histc(occinc,branges);
        inc_t(:,it,i)=h(1:10);
        %         histc(occinc,56000:16000:200000)
        
        avgrent_t(:,:,i)=avgrentstats{batchind(:,1)==i & batchind(:,2)==1};
        inc_wghavg(it,:)=sum(reshape(inc_t(:,it,:),10,4).*...
            repmat(branges(1:10)',1,4),1)./sum(reshape(inc_t(:,it,:),10,4),1);

    end
    ltmap=cat(2,LTmap{:,batchind(:,1)==i & batchind(:,2)==1});
    htype_run(:,:,i)=histc(ltmap,1:HT);
    vacrates(i,:)=vacrlts(batchind(:,1)==i & batchind(:,2)==1,:);
    
    % land sale price, location, and time
    landprice(:,i)=landsales(:,batchind(:,1)==i & batchind(:,2)==1);
    landtime(:,i)=landsaletime(:,batchind(:,1)==i & batchind(:,2)==1);
    [farmsale,ia,ic]=unique(landprice(:,i),'stable');
    farmsale=farmsale(farmsale~=0);
    saletime=landtime(ia,i);
    saletime=saletime(saletime~=0);
    [row,coastd]=ind2sub([NLENGTH NWIDTH],ia);
    landsale_run(i)=mat2cell([farmsale saletime (NWIDTH-coastd(2:length(coastd)))...
        dist2cbd(ia(2:length(ia)))],length(farmsale),4);
    
    % Rent sort
    for ht=1:HT
        rentmat(i,:,ht)=avgrent_t(ht,:,i);
    end
end


cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\figs\alt_stormclim
brlabel={'$40-$56','$56-$72','$72-$88','$88-$104','$104-$120','$120-$136',...
    '$13-$152','$152-$168','$168-$184','$184-$200'};
strmclim{1}='MidAt';
strmclim{2}='NC';
strmclim{3}='FL';
strmclim{4}='TX';
htlabel={'1','2','3','4','5','6','7','8'};

% hh3=figure;
% set(hh3,'color','white')

for i=1:4
    % income distribution
    hh1=figure;
    set(hh1,'color','white','Visible','off')
    set(hh1,'colormap',jet(10))
    area(TSTART+1:TMAX,(inc_t(:,TSTART+1:TMAX,i)./...
        repmat(max(sum(inc_t(:,TSTART+1:TMAX,i),1),1),10,1))')
    lcolorbar(brlabel,'TitleString','Income($ 10^3)');
    axis([TSTART+1 TMAX 0 1])
    title(strmclim(i))
    ylabel('Proportion Income Class')
    xlabel('Time Step')
    saveas(hh1,sprintf('Income_class_run1_%s',strmclim{i}),'jpg')

    % housing type
    hh4=figure;
    set(hh4,'color','white','Visible','off')
    set(hh4,'colormap',jet(8))
    area(TSTART+1:TMAX,(htype_run(:,TSTART+1:TMAX,i))')
%     area(TSTART+1:TMAX,(htype_run(:,TSTART+1:TMAX,i)./...
%         repmat(sum(htype_run(:,TSTART+1:TMAX,i),1),HT,1))')
    lcolorbar(htlabel,'ColorAlignment','center','TitleString','Housing Types')
    title(strmclim(i))
    ylabel('Housing Stock by Type')
    xlabel('Time Step')
    saveas(hh4,sprintf('Housing_Types_run1_%s',strmclim{i}),'jpg')
    
    % land sales space-time plot
    hh3=figure;
    set(hh3,'color','white','Visible','off')
    ls_data=landsale_run{i};
%     plot(ls_data(:,2)-TSTART,ls_data(:,3),'Marker','.','Color',reshape(rgb,length(rgb),3))
    ls_surf=zeros(length(TSTART+1:TMAX),NWIDTH);
    ls_surf(sub2ind(size(ls_surf),ls_data(:,2)-TSTART,ls_data(:,3)))=ls_data(:,1);
    imagesc(ls_surf)
    cmap=get(gcf,'colormap');
    cmap(1,:)=[1,1,1];
    set(hh3,'colormap',cmap)
    axis xy
%     set(gca,'XDir','reverse')
    title(sprintf('Land Sales Space-Time Plot, %s',strmclim{i}))
    xlabel('Distance from Coast')
    ylabel('Time of Land Sale')
    set(gca,'clim',[0 13000])
    colorbar
    saveas(hh3,sprintf('Land_Sales_run1_%s',strmclim{i}),'jpg')
end
% median rents
    hh2=figure;
    set(hh2,'color','white','Visible','off')
    plot(10:30,reshape(median(avgrent_t(:,10:30,:),1),21,4),'-')
    legend('Mid-Atl.','NC','FL','TX')
    ylabel('Median Housing Rent')
    xlabel('Time Step')
    saveas(hh2,sprintf('Median_rents_run1_%s',strmclim{i}),'jpg')


%% Create figures and results

cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\figs\alt_stormclim

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
strmclim{1}='MidAt';
strmclim{2}='NC';
strmclim{3}='FL';
strmclim{4}='TX';
h2_1=figure;
set(h2_1, 'Color','white','OuterPosition',[1,1,400,700],'Visible','off');

% plot housing types
for j=1:ERUNS
    h1=figure;
    set(h1, 'Color','white','Position',[1,1,700,700],'Visible','off');
    subplot(2,2,1);
    imagesc(AVGMAP(:,:,15,j));
    CMAP=colormap;
    CMAP(1,:)=[1 1 1];
    colormap(h1,CMAP);
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    title('t=15');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,2,2);
    imagesc(AVGMAP(:,:,20,j));
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    title('t=20');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,2,3);
    imagesc(AVGMAP(:,:,25,j));
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    title('t=25');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,2,4);
    imagesc(AVGMAP(:,:,30,j));
    set(gca,'xticklabel',[],'yticklabel',[])
%     set(gca,'Visible','off')
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
    CMAP=colormap;
    CMAP(1,:)=[1 1 1];
    colormap(h2,CMAP);
    set(gca,'clim',[40000 200000])
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    colorbar;
    title('Median Consumer Income');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,1,2);
    imagesc(reshape(varinc(:,j),NLENGTH,NWIDTH));
%     set(gca,'Visible','off')
    set(gca,'clim',[0 3600000000])
    set(gca,'xticklabel',[],'yticklabel',[])
    colorbar;
    title('Variance in Consumer Income');
    set(get(gca,'Title'),'Visible','on')
    saveas(h2,sprintf('final_consumer_income_%d',j),'jpg')
    
    figure(h2_1);
    subplot(ERUNS,1,j);
    hist(subincomemap(subincomemap(:,30)~=0,30),20);
    axis([40000 200000 0 200])
    if j < ERUNS
        set(gca,'xticklabel',[])
        title(sprintf('Income Distribution %s',strmclim{j}))
    else
        xlabel('Consumer Income')
        title(sprintf('Income Distribution %s',strmclim{j}))
    end
    
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
    CMAP=colormap;
    CMAP(1,:)=[1 1 1];
    colormap(h3,CMAP);
    set(gca,'clim',[9 TMAX])
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    colorbar;
    title('Average Build Time');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,1,2);
    imagesc(reshape(varbt(:,j),NLENGTH,NWIDTH));
%     set(gca,'Visible','off')
    set(gca,'clim',[0 TMAX])
    set(gca,'xticklabel',[],'yticklabel',[])
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
%     CMAP=colormap;
%     CMAP(1,:)=[1 1 1];
%     colormap(h5,CMAP);
    set(gca,'clim',[0 13000])
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    colorbar;
    title('Average Land Sales');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,1,2);
    imagesc(reshape(varlandsale(:,j),NLENGTH,NWIDTH));
%     set(gca,'Visible','off')
    set(gca,'clim',[0 12000000])
    set(gca,'xticklabel',[],'yticklabel',[])
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
    CMAP=colormap;
    CMAP(1,:)=[1 1 1];
    colormap(h6,CMAP);
    set(gca,'clim',[0 0.3])
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    colorbar;
    title('Mean Consumer Amenity Preference');
    set(get(gca,'Title'),'Visible','on')    
    subplot(2,1,2);
    imagesc(reshape(varampref(:,j),NLENGTH,NWIDTH));
%     set(gca,'Visible','off')
    set(gca,'clim',[0 0.03])
    set(gca,'xticklabel',[],'yticklabel',[])
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
    legend('HT1','HT2','HT3','HT4','HT5','HT6','HT7','HT8','Location','NorthWest');
    saveas(h7,sprintf('avgrents_%d',j),'jpg')
end

inc=cell(1,ERUNS);
incomedist=zeros(11,4);
htdist=zeros(HT,ERUNS);
rentdist=zeros(HT,ERUNS);
for q=1:ERUNS
    iruns=batchruns{q};
    subincomemap=incomemap(:,batchruns{q});
    inc(q)=mat2cell(subincomemap(subincomemap(:,30)~=0,30),...
        length(find(subincomemap(:,30)~=0)),1);
    incomedist(:,q)=hist(cat(1,inc{q}),[39999 56000 72000 88000 104000 ...
        120000 136000 152000 168000 184000 200001]);
    
    htdist(:,q)=mean(cat(2,numltrlts{30,iruns}),2);
    rentdist(:,q)=mean(avgrents(:,30,iruns),3);
    
end


%%

%%%% Summary plots
% Build Time
h_bt=figure;
set(h_bt,'Color','white','visible','off');
for k=1:ERUNS
    % xaxis = lvset
    %yaxis = am_slope
    subtightplot(5,5,k,0.01,0.1,0.05);
    imagesc(reshape(avgbt(:,k),NLENGTH,NWIDTH))
    set(gca,'xticklabel',[],'yticklabel',[])
    if isempty(find(ismember(1:5,k),1))==0
        title(sprintf('lv(%0.2f)',double(batchparms_unq(k,2))))
    end
    icol=ismember([1 6 11 16 21],k);
    if isempty(find(icol,1))==0
        ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
    end
    if k==23
        xlabel('Build Time (baseline: (0.1,1.0)','FontSize',14)
    end
end
saveas(h_bt,'summary_build_time','jpg')

h_ent=figure;
set(h_ent, 'Color','white','Visible','off');
for b=1:ERUNS
    subtightplot(5,5,b,0.01,0.1,0.05);
    plot(10:30,mean(htentropy(batchruns{b},10:30),1),'--k')
    hold on
    plot(10:30,min(htentropy(batchruns{b},10:30),[],1),'-b','LineWidth',0.5)
    plot(10:30,max(htentropy(batchruns{b},10:30),[],1),'-b','LineWidth',0.5)
    if b == 25
        set(gca,'yaxislocation','right')
    elseif b==23
        set(gca,'xticklabel',[],'yticklabel',[])
        xlabel('Avg. Entropy (baseline: (0.1,1.0)','FontSize',14)
    else
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if isempty(find(ismember(1:5,b),1))==0
        title(sprintf('lv(%0.2f)',double(batchparms_unq(b,2))))
    end
    icol=ismember([1 6 11 16 21],b);
    if isempty(find(icol,1))==0
        ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
    end 
end
saveas(h_ent,'summary_entropy','jpg')

% Average Map
h_avm=figure;
set(h_avm,'Color','white','visible','off');
for k=1:ERUNS
    % xaxis = lvset
    %yaxis = am_slope
    subtightplot(5,5,k,0.01,0.1,0.05);
    imagesc(AVGMAP(:,:,30,k))
    set(gca,'xticklabel',[],'yticklabel',[])
    if isempty(find(ismember(1:5,k),1))==0
        title(sprintf('lv(%0.2f)',double(batchparms_unq(k,2))))
    end
    icol=ismember([1 6 11 16 21],k);
    if isempty(find(icol,1))==0
        ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
    end
    if k==23
        xlabel('Avg. Map t=30 (baseline: (0.1,1.0)','FontSize',14)
    end
end
saveas(h_avm,'summary_avgmap','jpg')

% avg rents
avgrent_mat=cell(EXPTRUNS);
incdist_mat=cell(EXPTRUNS);
landsale_mat=cell(EXPTRUNS);
for jj=1:ERUNS
    % aggregate results for each parameter combo (n=30) and organize
    % according to plotting procedure to ensure the correct assignment of
    % results to parm settings.
    [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],jj);
    iruns=(matguide(:,1) == ir & matguide(:,2) == ic);
    avgrent_mat(ir,ic)=mat2cell(median(avgrents(:,TMAX,iruns),3),HT,1);
    incdist_mat(ir,ic)=mat2cell(histc(mdninc(mdninc(:,jj)~=0,jj),...
        linspace(20000,200000,20)),20,1);
    landsale_mat(ir,ic)=mat2cell(histc(unique(avglandsale(avglandsale(:,jj)~=0,jj)),...
        linspace(min(min(avglandsale)),max(max(avglandsale)),10)),10,1);
end
    
h_rent=figure;
set(h_rent,'Color','white','visible','off');
for a=1:ERUNS
    % plot avgrent histogram
    [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],a);
    subtightplot(5,5,a,0.01,0.1,0.05);
    bar(avgrent_mat{ir,ic})
    axis([0.5 8.5 0 22000]);
    if isempty(find(ismember([5 10 15 20],a),1))==0
         set(gca,'xticklabel',[],'yaxislocation','right')
    elseif a == 25
        set(gca,'yaxislocation','right')
    elseif a==23
        set(gca,'xticklabel',[],'yticklabel',[])
        xlabel('Avg. Rent t=30 (baseline: (0.1,1.0)','FontSize',14)
    else
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if isempty(find(ismember(1:5,a),1))==0
        title(sprintf('lv(%0.2f)',double(batchparms_unq(a,2))))
    end
    icol=ismember([1 6 11 16 21],a);
    if isempty(find(icol,1))==0
        ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
    end 
end

h_inc=figure;
set(h_inc,'visible','off');
for a=1:ERUNS
    % plot inc. dist. histogram
    [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],a);
    subtightplot(5,5,a,0.01,0.1,0.05);
    bar(linspace(20000,200000,20),incdist_mat{ir,ic},'histc')
    axis([10000 210000 0 750])
    if isempty(find(ismember([5 10 15 20],a),1))==0
         set(gca,'xticklabel',[],'yaxislocation','right')
    elseif a == 25
        set(gca,'yaxislocation','right')
    elseif a==23
        set(gca,'xticklabel',[],'yticklabel',[])
        xlabel('Income t=30 (baseline: (0.1,1.0)','FontSize',14)
    else
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if isempty(find(ismember(1:5,a),1))==0
        title(sprintf('lv(%0.2f)',double(batchparms_unq(a,2))))
    end
    icol=ismember([1 6 11 16 21],a);
    if isempty(find(icol,1))==0
        ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
    end 
end
saveas(h_inc,'summary_incdist','jpg')

h_ls=figure;
set(h_ls,'visible','off');
for a=1:ERUNS
    % plot inc. dist. histogram
    [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],a);
    subtightplot(5,5,a,0.01,0.1,0.05);
    bar(linspace(min(min(avglandsale)),max(max(avglandsale)),10),landsale_mat{ir,ic},'histc')
    axis([500 27000 0 40])
    if a == 25
        set(gca,'yaxislocation','right')
    elseif a==23
        set(gca,'xticklabel',[],'yticklabel',[])
        xlabel('Land Sales t=30 (baseline: (0.1,1.0)','FontSize',14)
    else
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if isempty(find(ismember(1:5,a),1))==0
        title(sprintf('lv(%0.2f)',double(batchparms_unq(a,2))))
    end
    icol=ismember([1 6 11 16 21],a);
    if isempty(find(icol,1))==0
        ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
    end 
end
saveas(h_ls,'summary_landsales','jpg')
