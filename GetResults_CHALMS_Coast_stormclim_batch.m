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
z = 1000*[0.00025    1.5000
    0.00025    2.5000
    0.0005    1.5000
    0.0005    2.5000
    0.0010    1.5000
    0.0010    2.5000
    0.0020    1.5000
    0.0020    2.5000];
% number of model iterations
MRUNS=30;
% number of experimental parameters
NPARMS=1;
% number of experimental variatons on NPARMS
EXPTRUNS=3;
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
cd X:\model_results\CHALMS_coast_032315
fnames=dir;
fnamescell=struct2cell(fnames);
% (2) change the prefix of the results file names
h=strncmp('coast_baseline',fnamescell(1,:),14);
hind=find(h==1);
% add precalculated distance matrix
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files
load DIST2CBD_east
cd X:\model_results\CHALMS_coast_032315
distpt=ceil(min(min(dist2cbd)):max(max(dist2cbd)));
density=zeros(NLENGTH*NWIDTH,length(hind));

pctdev=zeros(length(hind),TMAX);
htprob=cell(TMAX,length(hind));
htentropy=zeros(length(hind),TMAX);

htprob_cross=cell(TMAX,length(hind));
htentropy_cross=zeros(length(hind),TMAX);

LTmap=cell(TMAX,length(hind));
avgrentstats=cell(1,length(hind));
avgreturnstats=cell(1,length(hind));
avgexptreturnstats=cell(1,length(hind));
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
for mr=1:length(hind)   % MRUNS*EXPTRUNS
    h=strcmp(sprintf('coast_baseline%d_%d.mat',batchind(mr,1),...
        batchind(mr,2)),fnamescell(1,:));
    filename=fnamescell{1,h};
    load(filename)
    
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
%%
% % Calculate housing stock evenness and build average development map
for ie=1:EXPTRUNS
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
            ivarthresh=(devprob >= subavgthresh(iter));
            locmat=find(ivarthresh==1);
            testpctdev(ie,tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
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
    end
end
%% Create struct files of saved results
% spatial data at t=TMAX, size=[NLENGTH*NWIDTH,MRUNS*EXPTRUNS]
mapdata_store=struct('lot_types',{LTmap},'income_map',{incomemap},...
    'amenity_prefs',{amprefmap},'house_prices',{housepricemap},...
    'build_time',{btmap},'land_sales_price',{landsales},'land_sales_time',...
    {landsaletime},'average_maps',{AVGMAP});

% results by housing type over time, size=[HT,TMAX,MRUNS*EXPTRUNS]
numlts_store=cell(1,MRUNS*EXPTRUNS);
for ir=1:MRUNS*EXPTRUNS
    numlts_store(ir)=mat2cell(cat(2,numltrlts{:,ir}),HT,TMAX);
end
htdata_store=struct('num_lots',{numlts_store},'avg_rents',{avgrentstats},...
    'avg_return',{avgreturnstats},'avg_expt_ret',{avgexptreturnstats},...
    'house_stock_evenness',{htentropy});

% aggregate data across housing types
% build consumer income distribution
incdist_cell=cell(1,length(hind));
for j=1:length(hind)
    inc1=incomemap(:,j);
    bt1=btmap(:,j);
    inc_t=zeros(length(branges)-1,length(TSTART:TMAX));
    for it=TSTART:TMAX
        ibt=(bt1==it);
        occinc=unique(inc1(ibt));
        occinc=occinc(occinc~=0);
        if isempty(find(occinc,1))==1
            inc_t(:,it-9)=zeros(length(inc_t(:,1)),1);
            continue
        end
        h=histc(occinc,branges);
        inc_t(:,it-9)=h(1:length(branges)-1);
    end
    incdist_cell(j)=mat2cell(inc_t,length(branges)-1,length(TSTART:TMAX));
end
aggdata_store=struct('vacany_rates',{vacrlts},'income_distribution',{incdist_cell});

run_indices=struct('batch_indices',{batchind});

% save results_stormclim_amenityslope_batch_struct mapdata_store htdata_store aggdata_store run_indices

%%
landsale_run=cell(MRUNS,EXPTRUNS);
landsalerlts_time=zeros(TMAX,EXPTRUNS,7);
landsalerlts_coast=zeros(8,EXPTRUNS,7);
landsalerlts_cbd=zeros(11,EXPTRUNS,7);
saletime=cell(1,4);
incomedist=zeros(length(branges)-1,TMAX,EXPTRUNS);
hdistpct=zeros(length(branges)-1,TMAX,MRUNS);
plotpairs=cell(1,EXPTRUNS);
for i=1:EXPTRUNS
    iruns=batchruns{i};
    sublandsales=landsales(:,batchruns{i});
    sublandtime=landsaletime(:,batchruns{i});
    subincmap=incomemap(:,iruns);
    subbtmap=btmap(:,iruns);
    for ii=1:MRUNS
        for it=TSTART:TMAX
            ibt=(subbtmap(:,ii)==it);
            if isempty(find(ibt,1))==1
                hdistpct(:,it,ii)=hdistpct(:,it-1,ii);
                continue
            end
            hdist=histc(subincmap(ibt,ii),branges);
            if isempty(find(hdist,1))==1
                hdistpct(:,it,ii)=hdistpct(:,it-1,ii);
                continue
            end
            hdistpct(:,it,ii)=hdist(1:length(branges)-1)./sum(hdist(1:length(branges)-1));
        end

        % land sale price, location, and time
        [farmsale,ia,ic]=unique(sublandsales(:,ii),'stable');
        ia=ia(farmsale~=0);
        farmsale=farmsale(farmsale~=0);
        saletime=sublandtime(ia,ii);
        [row,coastd]=ind2sub([NLENGTH NWIDTH],ia);
        landsale_run(ii,i)=mat2cell([farmsale saletime (NWIDTH-coastd)...
            dist2cbd(ia)],length(farmsale),4);
    end
    incomedist(:,:,i)=mean(hdistpct,3);
    holdpairs=cat(1,landsale_run{:,i});
    plotpairs(i)=mat2cell(holdpairs(:,1:2),length(holdpairs),2);
    landsale_all=cat(1,landsale_run{:,i});
    for it=TSTART:TMAX
        isale=(landsale_all(:,2)==it);
        if isempty(find(isale==1,1))==1
            continue
        end
        [lsmu,lssigma,lsmuci,lssigmaci]=normfit(landsale_all(isale,1));
        lsse=lssigma/sqrt(length(landsale_all(isale,1)));
        landsalerlts_time(it,i,1)=lsmu;
        landsalerlts_time(it,i,2)=max(landsale_all(isale,1));
        landsalerlts_time(it,i,3)=min(landsale_all(isale,1));
        landsalerlts_time(it,i,4)=lssigma;
        landsalerlts_time(it,i,5:6)=lsmuci;
        landsalerlts_time(it,i,7)=lsse;
    end
    coastdistpt=unique(landsale_all(:,3));
    for ic=1:length(coastdistpt)
        isale=(landsale_all(:,3)==coastdistpt(ic));
        if isempty(find(isale==1,1))==1
            continue
        end
        [lsmu,lssigma,lsmuci,lssigmaci]=normfit(landsale_all(isale,1));
        lsse=lssigma/sqrt(length(landsale_all(isale,1)));
        landsalerlts_coast(ic,i,1)=lsmu;
        landsalerlts_coast(ic,i,2)=max(landsale_all(isale,1));
        landsalerlts_coast(ic,i,3)=min(landsale_all(isale,1));
        landsalerlts_coast(ic,i,4)=lssigma;
        landsalerlts_coast(ic,i,5:6)=lsmuci;
        landsalerlts_coast(ic,i,7)=lsse;
    end
    distpts=reshape(distpt(1:110),10,11);
    cbddist=zeros(2,11);
    cbddist(1,:)=distpts(1,:);
    cbddist(2,:)=distpts(10,:);
    for id=1:length(cbddist(1,:))
        isale=(landsale_all(:,4)>=cbddist(1,id) & landsale_all(:,4)<=cbddist(2,id));
        if isempty(find(isale==1,1))==1
            continue
        end
        [lsmu,lssigma,lsmuci,lssigmaci]=normfit(landsale_all(isale,1));
        lsse=lssigma/sqrt(length(landsale_all(isale,1)));
        landsalerlts_cbd(id,i,1)=lsmu;
        landsalerlts_cbd(id,i,2)=max(landsale_all(isale,1));
        landsalerlts_cbd(id,i,3)=min(landsale_all(isale,1));
        landsalerlts_cbd(id,i,4)=lssigma;
        landsalerlts_cbd(id,i,5:6)=lsmuci;
        landsalerlts_cbd(id,i,7)=lsse;
    end
end
land_sales=struct('landsales_time',{landsalerlts_time},'landsales_coast',...
    {landsalerlts_coast},'landsales_cbd',{landsalerlts_cbd},'landsales_all',{landsale_run});
save results_coast_baseline_batch_struct mapdata_store htdata_store aggdata_store run_indices land_sales


% %%% Plot land price trajectories
% landsale_traj=landsalerlts_time(:,:,1);
% landsale_traj=landsale_traj(12:30,:)';
% landsale_traj_025=landsale_traj([1 8 15 22],:);
% landsale_traj_050=landsale_traj([2 9 16 23],:);
% landsale_traj_075=landsale_traj([3 10 17 24],:);
% landsale_traj_100=landsale_traj([4 11 18 25],:);
% landsale_traj_125=landsale_traj([5 12 19 26],:);
% landsale_traj_150=landsale_traj([6 13 20 27],:);
% landsale_traj_175=landsale_traj([7 14 21 28],:);
% 
% landsale_traj_patches=zeros(7,2*length(landsale_traj(1,:))+1);
% patchx=[12:30 30:-1:12 12];
% landsale_traj_patches(1,:)=[max(landsale_traj_025,[],1) min(landsale_traj_025,[],1) max(landsale_traj_025(:,1))];
% landsale_traj_patches(2,:)=[max(landsale_traj_050,[],1) min(landsale_traj_050,[],1) max(landsale_traj_050(:,1))];
% landsale_traj_patches(3,:)=[max(landsale_traj_075,[],1) min(landsale_traj_075,[],1) max(landsale_traj_075(:,1))];
% landsale_traj_patches(4,:)=[max(landsale_traj_100,[],1) min(landsale_traj_100,[],1) max(landsale_traj_100(:,1))];
% landsale_traj_patches(5,:)=[max(landsale_traj_125,[],1) min(landsale_traj_125,[],1) max(landsale_traj_125(:,1))];
% landsale_traj_patches(6,:)=[max(landsale_traj_150,[],1) min(landsale_traj_150,[],1) max(landsale_traj_150(:,1))];
% landsale_traj_patches(7,:)=[max(landsale_traj_175,[],1) min(landsale_traj_175,[],1) max(landsale_traj_175(:,1))];
% 
% figure
% patch(patchx,landsale_traj_patches(7,:),'r')
% hold on
% patch(patchx,landsale_traj_patches(6,:),'y')
% patch(patchx,landsale_traj_patches(5,:),'g')
% patch(patchx,landsale_traj_patches(4,:),'k')
% patch(patchx,landsale_traj_patches(3,:),'b')
% patch(patchx,landsale_traj_patches(2,:),'c')
% patch(patchx,landsale_traj_patches(1,:),'m')
% 
% 
% 
% 
% 
% 
% plot(12:30,min(landsale_traj_025,[],1),'-m','LineWidth',1.5)
% hold on
% plot(12:30,max(landsale_traj_025,[],1),'-m','LineWidth',1.5)
% plot(12:30,min(landsale_traj_050,[],1),'-c','LineWidth',1.5)
% plot(12:30,max(landsale_traj_050,[],1),'-c','LineWidth',1.5)
% plot(12:30,min(landsale_traj_075,[],1),'-b','LineWidth',1.5)
% plot(12:30,max(landsale_traj_075,[],1),'-b','LineWidth',1.5)
% plot(12:30,min(landsale_traj_100,[],1),'-k','LineWidth',1.5)
% plot(12:30,max(landsale_traj_100,[],1),'-k','LineWidth',1.5)
% plot(12:30,min(landsale_traj_125,[],1),'-g','LineWidth',1.5)
% plot(12:30,max(landsale_traj_125,[],1),'-g','LineWidth',1.5)
% plot(12:30,min(landsale_traj_150,[],1),'-y','LineWidth',1.5)
% plot(12:30,max(landsale_traj_150,[],1),'-y','LineWidth',1.5)
% plot(12:30,min(landsale_traj_175,[],1),'-r','LineWidth',1.5)
% plot(12:30,max(landsale_traj_175,[],1),'-r','LineWidth',1.5)
% legend('0.025','0.025','0.05','0.05','0.075','0.075','0.10','0.10','0.125',...
%     '0.125','0.150','0.150','0.175','0.175','Location','northwest')
%%
runnamelabel={'MidAtl','NC','FL','TX'};
tset=mat2cell(TSTART:TMAX,1,length(TSTART:TMAX));
cdistlabel=mat2cell(((coastdistpt.*cell2mile)-coastdistpt(1)*cell2mile)',1,8);
cbddistlabel=mat2cell(cbddist(1,:),1,11);
rentlabel={'Avg Rents'};
retlabel={'Avg Returns'};
htutillabel={'Avg Utility'};
novaclabel={'No Vacancies'};
htincomelabel={'Income per ht'};
bidlabel={'Share of Bids'};
statlabel={'Low 95 CI','Up 95 CI','SE'}';
lotlabel={'Lots Per Type'};
lotsumlabel={'Total Lots'};
lotsizelabel={'Avg Lot Size'};
hpaclabel={'Houses Per Acre'};
vaclabel={'Vacancy Rates'};
utillabel={'Pct of Max Utility'};
displabel={'Mean Dispersal'};
devlabel={'Pct Developed Area'};
biglotslabel={'Pct Development >= 1 acre'};
zonedenlabel={'Avg. Lot Size '};
distlabel={'Distance from CBD'};
distlabel_c={'Distance from Coast'};
incomelabel={'Median Income of Occupants'};
incdistlabel={'Median Income Distribution (% of consumers)'};
northincomelabel={'Mean Income of Occupants in Unzoned Region'};
southincomelabel={'Mean Income of Occupants in Zoned Region'};
coastincomelabel={'Mean Income of Occupants in Coastal Zone'};
middleincomelabel={'Mean Income of Occupants in Midland Zone'};
inlandincomelabel={'Mean Income of Occupants in Inland Zone'};
outincomelabel={'Mean Income of Ex-occupants'};
zonedpctdevlabel={'Pct. Dev. in Zoned Region'};
unzonedpctdevlabel={'Pct. Dev. in Unzoned Region'};
coastalpctdevlabel={'Pct. Dev. in Coastal Region'};
middlepctdevlabel={'Pct. Dev. in Midland Region'};
inlandpctdevlabel={'Pct. Dev. in Inland Region'};
maxdistlabel={'Max Linear Distance of Development'};
zonepricelabel={'Mean Land Price by Zoning Region'};
coastpricelabel={'Mean Land Price by Coast/Midland/Inland Region'};
northlabel={'North'};
southlabel={'South'};
coastlabel={'Coastal'};
middlelabel={'Midland'};
inlandlabel={'Inland'};
% psoldlabel={'Mean Agricultural Return'};
pfarmlabel={'Mean Agricultural Return'};
farmsoldlabel={'Sold'};
farmunsoldlabel={'Unsold'};
timelabel={'Time Step'};
lottypelabel={'Lot Type'};
farmeridlabel={'Farmer ID'};
selltimelabel={'Sell Time'};
selpricelabel={'Land Price'};
lottypeset=(1:(HT))';
farminfolabel={'Farm Info'};
planddistlabel={'Distance from CBD (mi) vs. Land Price ($/ac)'};
planddistclabel={'Distance from Coast (mi) vs. Land Price ($/ac)'};
plandtimelabel={'Time Step of Land Sale vs. Land Price ($/ac)'};
tlabel={'Time Step'};
dlabel={'Distance'};
pllabel={'Land Price'};
plandstatslabel={'Land Price Stats'};
meanlabel={'Mean'};
maxlabel={'Max'};
minlabel={'Min'};
sigmalabel={'Sigma'};
offerlabel={'Avg House Offering'};
ideallabel={'Utility Max House Offering'};
proflabel={'Profit Max House Offering'};
difflabel={'Difference'};
dscptlabel={'Descriptive Stats'};
carrycostlabel={'Carrying Cost'};
binlabel={'Bin Min. ($)'};
landvaluelabel={'Initial Land Value Setting'};

% %% Write time step file %%
% cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\results\amenityslope_landvalues
% resultsfile=('Results_CHALMS_amenityslope_landvalues.xlsx');
% for j=1:EXPTRUNS
%     xlswrite(resultsfile,timelabel,sprintf('Time Step Stats %s',runnamelabel{j}),'C1');
%     xlswrite(resultsfile,tset{:},sprintf('Time Step Stats %s',runnamelabel{j}),'C2');
%     xlswrite(resultsfile,lottypelabel,sprintf('Time Step Stats %s',runnamelabel{j}),'B2');
%     xlswrite(resultsfile,lottypeset,sprintf('Time Step Stats %s',runnamelabel{j}),'B3');
%     xlswrite(resultsfile,rentlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A3');
%     rentdump=reshape(cell2mat(avgrentstats(batchruns{j})),HT,TMAX,MRUNS);
%     xlswrite(resultsfile,mean(rentdump(:,TSTART:TMAX,:),3),sprintf('Time Step Stats %s',runnamelabel{j}),'C3');
%     
%     xlswrite(resultsfile,lotlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A12');
%     xlswrite(resultsfile,lottypeset,sprintf('Time Step Stats %s',runnamelabel{j}),'B12');
%     ltdump=reshape(cell2mat(numltrlts(:,batchruns{j})),HT,TMAX,MRUNS);
%     xlswrite(resultsfile,mean(ltdump(:,TSTART:TMAX,:),3),sprintf('Time Step Stats %s',runnamelabel{j}),'C12');
%     xlswrite(resultsfile,lotsumlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'B20');
%     xlswrite(resultsfile,sum(mean(ltdump(:,TSTART:TMAX,:),3),1),sprintf('Time Step Stats %s',runnamelabel{j}),'C20');
%     
%     xlswrite(resultsfile,dscptlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A22');
%     xlswrite(resultsfile,timelabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A23');
%     xlswrite(resultsfile,tset{:},sprintf('Time Step Stats %s',runnamelabel{j}),'B23');
%     xlswrite(resultsfile,devlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A24');
%     xlswrite(resultsfile,mean(pctdev(batchruns{j},TSTART:TMAX),1),sprintf('Time Step Stats %s',runnamelabel{j}),'B24');
%     xlswrite(resultsfile,vaclabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A25');
%     xlswrite(resultsfile,avgvac(j,:),sprintf('Time Step Stats %s',runnamelabel{j}),'B25');
%     xlswrite(resultsfile,incdistlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A27');
%     xlswrite(resultsfile,binlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A28');
%     xlswrite(resultsfile,branges',sprintf('Time Step Stats %s',runnamelabel{j}),'B28');
%     xlswrite(resultsfile,incomedist(:,TSTART:TMAX,j),sprintf('Time Step Stats %s',runnamelabel{j}),'C28');
% 
%     xlswrite(resultsfile,plandstatslabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A39');
%     xlswrite(resultsfile,plandtimelabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A40');
%     xlswrite(resultsfile,tlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A41');
%     xlswrite(resultsfile,meanlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A42');
%     xlswrite(resultsfile,maxlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A43');
%     xlswrite(resultsfile,minlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A44');
%     xlswrite(resultsfile,sigmalabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A45');
%     xlswrite(resultsfile,statlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A46');
%     xlswrite(resultsfile,tset{:},sprintf('Time Step Stats %s',runnamelabel{j}),'B41');
%     xlswrite(resultsfile,reshape(landsalerlts_time(TSTART:TMAX,j,:),length(TSTART:TMAX),7)',...
%         sprintf('Time Step Stats %s',runnamelabel{j}),'B42');
%     
%     xlswrite(resultsfile,plandstatslabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A50');
%     xlswrite(resultsfile,planddistclabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A51');
%     xlswrite(resultsfile,dlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A52');
%     xlswrite(resultsfile,meanlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A53');
%     xlswrite(resultsfile,maxlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A54');
%     xlswrite(resultsfile,minlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A55');
%     xlswrite(resultsfile,sigmalabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A56');
%     xlswrite(resultsfile,statlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A57');
%     xlswrite(resultsfile,cdistlabel{:},sprintf('Time Step Stats %s',runnamelabel{j}),'B52');
%     xlswrite(resultsfile,reshape(landsalerlts_coast(:,j,:),8,7)',...
%         sprintf('Time Step Stats %s',runnamelabel{j}),'B53');
%     
%     xlswrite(resultsfile,plandstatslabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A61');
%     xlswrite(resultsfile,planddistclabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A62');
%     xlswrite(resultsfile,dlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A63');
%     xlswrite(resultsfile,meanlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A64');
%     xlswrite(resultsfile,maxlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A65');
%     xlswrite(resultsfile,minlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A66');
%     xlswrite(resultsfile,sigmalabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A67');
%     xlswrite(resultsfile,statlabel,sprintf('Time Step Stats %s',runnamelabel{j}),'A68');
%     xlswrite(resultsfile,cbddistlabel{:},sprintf('Time Step Stats %s',runnamelabel{j}),'B63');
%     xlswrite(resultsfile,reshape(landsalerlts_cbd(:,j,:),11,7)',...
%         sprintf('Time Step Stats %s',runnamelabel{j}),'B64');
% end

%% Single run comparisons
% irun1=(batchind(:,2)==1);
% ls1=landsales(:,batchind(:,2)==1);
% bt1=btmap(:,batchind(:,2)==1);
% inc1=incomemap(:,batchind(:,2)==1);
% 
% inc_t=zeros(10,TMAX,4);
% avgrent_t=zeros(HT,TMAX,4);
% vacrates=zeros(4,21);
% landprice=zeros(NCELLS,4);
% landtime=zeros(NCELLS,4);
% landsale_run=cell(1,4);
% saletime=cell(1,4);
% branges=[40000:16000:184000 200001];
% htype_run=zeros(HT,TMAX,4);
% rentmat=zeros(4,TMAX,HT);
% inc_wghavg=zeros(TMAX,4);
% for i=1:4
%     for it=TSTART:TMAX
%         ibt=(bt1(:,i)==it);
%         occinc=unique(inc1(ibt,i));
%         occinc=occinc(occinc~=0);
%         h=histc(occinc,branges);
%         inc_t(:,it,i)=h(1:10);
%         %         histc(occinc,56000:16000:200000)
%         
%         avgrent_t(:,:,i)=avgrentstats{batchind(:,1)==i & batchind(:,2)==1};
%         inc_wghavg(it,:)=sum(reshape(inc_t(:,it,:),10,4).*...
%             repmat(branges(1:10)',1,4),1)./sum(reshape(inc_t(:,it,:),10,4),1);
% 
%     end
%     ltmap=cat(2,LTmap{:,batchind(:,1)==i & batchind(:,2)==1});
%     htype_run(:,:,i)=histc(ltmap,1:HT);
%     vacrates(i,:)=vacrlts(batchind(:,1)==i & batchind(:,2)==1,:);
%     
%     % land sale price, location, and time
%     landprice(:,i)=landsales(:,batchind(:,1)==i & batchind(:,2)==1);
%     landtime(:,i)=landsaletime(:,batchind(:,1)==i & batchind(:,2)==1);
%     [farmsale,ia,ic]=unique(landprice(:,i),'stable');
%     farmsale=farmsale(farmsale~=0);
%     saletime=landtime(ia,i);
%     saletime=saletime(saletime~=0);
%     [row,coastd]=ind2sub([NLENGTH NWIDTH],ia);
%     landsale_run(i)=mat2cell([farmsale saletime (NWIDTH-coastd(2:length(coastd)))...
%         dist2cbd(ia(2:length(ia)))],length(farmsale),4);
%     
%     % Rent sort
%     for ht=1:HT
%         rentmat(i,:,ht)=avgrent_t(ht,:,i);
%     end
% end
% 
% 
% cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\figs\alt_stormclim
% brlabel={'$40-$56','$56-$72','$72-$88','$88-$104','$104-$120','$120-$136',...
%     '$13-$152','$152-$168','$168-$184','$184-$200'};
% strmclim{1}='MidAt';
% strmclim{2}='NC';
% strmclim{3}='FL';
% strmclim{4}='TX';
% htlabel={'1','2','3','4','5','6','7','8'};
% 
% % hh3=figure;
% % set(hh3,'color','white')
% 
% for i=1:4
%     % income distribution
%     hh1=figure;
%     set(hh1,'color','white','Visible','off')
%     set(hh1,'colormap',jet(10))
%     area(TSTART+1:TMAX,(inc_t(:,TSTART+1:TMAX,i)./...
%         repmat(max(sum(inc_t(:,TSTART+1:TMAX,i),1),1),10,1))')
%     lcolorbar(brlabel,'TitleString','Income($ 10^3)');
%     axis([TSTART+1 TMAX 0 1])
%     title(strmclim(i))
%     ylabel('Proportion Income Class')
%     xlabel('Time Step')
%     saveas(hh1,sprintf('Income_class_run1_%s',strmclim{i}),'jpg')
% 
%     % housing type
%     hh4=figure;
%     set(hh4,'color','white','Visible','off')
%     set(hh4,'colormap',jet(8))
%     area(TSTART+1:TMAX,(htype_run(:,TSTART+1:TMAX,i))')
% %     area(TSTART+1:TMAX,(htype_run(:,TSTART+1:TMAX,i)./...
% %         repmat(sum(htype_run(:,TSTART+1:TMAX,i),1),HT,1))')
%     lcolorbar(htlabel,'ColorAlignment','center','TitleString','Housing Types')
%     title(strmclim(i))
%     ylabel('Housing Stock by Type')
%     xlabel('Time Step')
%     saveas(hh4,sprintf('Housing_Types_run1_%s',strmclim{i}),'jpg')
%     
%     % land sales space-time plot
%     hh3=figure;
%     set(hh3,'color','white','Visible','off')
%     ls_data=landsale_run{i};
% %     plot(ls_data(:,2)-TSTART,ls_data(:,3),'Marker','.','Color',reshape(rgb,length(rgb),3))
%     ls_surf=zeros(length(TSTART+1:TMAX),NWIDTH);
%     ls_surf(sub2ind(size(ls_surf),ls_data(:,2)-TSTART,ls_data(:,3)))=ls_data(:,1);
%     imagesc(ls_surf)
%     cmap=get(gcf,'colormap');
%     cmap(1,:)=[1,1,1];
%     set(hh3,'colormap',cmap)
%     axis xy
% %     set(gca,'XDir','reverse')
%     title(sprintf('Land Sales Space-Time Plot, %s',strmclim{i}))
%     xlabel('Distance from Coast')
%     ylabel('Time of Land Sale')
%     set(gca,'clim',[0 13000])
%     colorbar
%     saveas(hh3,sprintf('Land_Sales_run1_%s',strmclim{i}),'jpg')
% end
% % median rents
%     hh2=figure;
%     set(hh2,'color','white','Visible','off')
%     plot(10:30,reshape(median(avgrent_t(:,10:30,:),1),21,4),'-')
%     legend('Mid-Atl.','NC','FL','TX')
%     ylabel('Median Housing Rent')
%     xlabel('Time Step')
%     saveas(hh2,sprintf('Median_rents_run1_%s',strmclim{i}),'jpg')
% 
% 
%% Create figures and results

cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\figs\coast_baseline

% resultsfile='baseline_results_070314.txt';
% save(resultsfile,'avgvac','avgpctdev','-ascii')
mdninc=zeros(NLENGTH*NWIDTH,EXPTRUNS);
varinc=zeros(NLENGTH*NWIDTH,EXPTRUNS);
avgbt=zeros(NLENGTH*NWIDTH,EXPTRUNS);
varbt=zeros(NLENGTH*NWIDTH,EXPTRUNS);
avglandsale=zeros(NLENGTH*NWIDTH,EXPTRUNS);
varlandsale=zeros(NLENGTH*NWIDTH,EXPTRUNS);
avgampref=zeros(NLENGTH*NWIDTH,EXPTRUNS);
varampref=zeros(NLENGTH*NWIDTH,EXPTRUNS);
strmclim{1}='MidAt';
strmclim{2}='NC';
strmclim{3}='FL';
strmclim{4}='TX';
h2_1=figure;
set(h2_1, 'Color','white','OuterPosition',[1,1,400,700],'Visible','off');

% plot housing types
for j=1:EXPTRUNS
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
    
%     figure(h2_1);
%     orient tall
%     set(h2_1, 'Color','white','OuterPosition',[1,1,400,700],'Visible','off');
%     subplot(EXPTRUNS,1,j);
% %     hist(subincomemap(subincomemap(:,30)~=0,30),20);
%     hist(mdninc(isubdev,j),20);
%     axis([40000 200000 0 700])
%     if j < EXPTRUNS
%         set(gca,'xticklabel',[])
%         title(sprintf('Income Distribution %s',strmclim{j}))
%     else
%         xlabel('Consumer Income')
%         title(sprintf('Income Distribution %s',strmclim{j}))
%     end
%     saveas(h2_1,'consumer_income_distributions','jpg')
    
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
    set(gca,'clim',[0 20000])
%     set(gca,'Visible','off')
    set(gca,'xticklabel',[],'yticklabel',[])
    colorbar;
    title('Average Land Sales');
    set(get(gca,'Title'),'Visible','on')
    subplot(2,1,2);
    imagesc(reshape(varlandsale(:,j),NLENGTH,NWIDTH));
%     set(gca,'Visible','off')
    set(gca,'clim',[0 20000000])
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

% inc=cell(1,ERUNS);
% incomedist=zeros(11,4);
% htdist=zeros(HT,ERUNS);
% rentdist=zeros(HT,ERUNS);
% for q=1:EXPTRUNS
%     iruns=batchruns{q};
%     subincomemap=incomemap(:,batchruns{q});
% %     inc(q)=mat2cell(subincomemap(subincomemap(:,30)~=0,30),...
% %         length(find(subincomemap(:,30)~=0)),1);
%     inc(q)=mat2cell(mdninc(mdninc(:,q)~=0,q),length(find(mdninc(:,q)~=0)),1);
%     incomedist(:,q)=hist(cat(1,inc{q}),[39999 56000 72000 88000 104000 ...
%         120000 136000 152000 168000 184000 200001]);
%     
%     htdist(:,q)=mean(cat(2,numltrlts{30,iruns}),2);
%     rentdist(:,q)=mean(avgrents(:,30,iruns),3);
%     
% end


% %% Plots for multiple experimental parameters
% 
% %%%% Summary plots
% % Build Time
% h_bt=figure;
% set(h_bt,'Color','white','visible','off');
% for k=1:ERUNS
%     % xaxis = lvset
%     %yaxis = am_slope
%     subtightplot(5,5,k,0.01,0.1,0.05);
%     imagesc(reshape(avgbt(:,k),NLENGTH,NWIDTH))
%     set(gca,'xticklabel',[],'yticklabel',[])
%     if isempty(find(ismember(1:5,k),1))==0
%         title(sprintf('lv(%0.2f)',double(batchparms_unq(k,2))))
%     end
%     icol=ismember([1 6 11 16 21],k);
%     if isempty(find(icol,1))==0
%         ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
%     end
%     if k==23
%         xlabel('Build Time (baseline: (0.1,1.0)','FontSize',14)
%     end
% end
% saveas(h_bt,'summary_build_time','jpg')
% 
% h_ent=figure;
% set(h_ent, 'Color','white','Visible','off');
% for b=1:ERUNS
%     subtightplot(5,5,b,0.01,0.1,0.05);
%     plot(10:30,mean(htentropy(batchruns{b},10:30),1),'--k')
%     hold on
%     plot(10:30,min(htentropy(batchruns{b},10:30),[],1),'-b','LineWidth',0.5)
%     plot(10:30,max(htentropy(batchruns{b},10:30),[],1),'-b','LineWidth',0.5)
%     if b == 25
%         set(gca,'yaxislocation','right')
%     elseif b==23
%         set(gca,'xticklabel',[],'yticklabel',[])
%         xlabel('Avg. Entropy (baseline: (0.1,1.0)','FontSize',14)
%     else
%         set(gca,'xticklabel',[],'yticklabel',[])
%     end
%     if isempty(find(ismember(1:5,b),1))==0
%         title(sprintf('lv(%0.2f)',double(batchparms_unq(b,2))))
%     end
%     icol=ismember([1 6 11 16 21],b);
%     if isempty(find(icol,1))==0
%         ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
%     end 
% end
% saveas(h_ent,'summary_entropy','jpg')
% 
% % Average Map
% h_avm=figure;
% set(h_avm,'Color','white','visible','off');
% for k=1:ERUNS
%     % xaxis = lvset
%     %yaxis = am_slope
%     subtightplot(5,5,k,0.01,0.1,0.05);
%     imagesc(AVGMAP(:,:,30,k))
%     set(gca,'xticklabel',[],'yticklabel',[])
%     if isempty(find(ismember(1:5,k),1))==0
%         title(sprintf('lv(%0.2f)',double(batchparms_unq(k,2))))
%     end
%     icol=ismember([1 6 11 16 21],k);
%     if isempty(find(icol,1))==0
%         ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
%     end
%     if k==23
%         xlabel('Avg. Map t=30 (baseline: (0.1,1.0)','FontSize',14)
%     end
% end
% saveas(h_avm,'summary_avgmap','jpg')
% 
% % avg rents
% avgrent_mat=cell(EXPTRUNS);
% incdist_mat=cell(EXPTRUNS);
% landsale_mat=cell(EXPTRUNS);
% for jj=1:ERUNS
%     % aggregate results for each parameter combo (n=30) and organize
%     % according to plotting procedure to ensure the correct assignment of
%     % results to parm settings.
%     [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],jj);
%     iruns=(matguide(:,1) == ir & matguide(:,2) == ic);
%     avgrent_mat(ir,ic)=mat2cell(median(avgrents(:,TMAX,iruns),3),HT,1);
%     incdist_mat(ir,ic)=mat2cell(histc(mdninc(mdninc(:,jj)~=0,jj),...
%         linspace(20000,200000,20)),20,1);
%     landsale_mat(ir,ic)=mat2cell(histc(unique(avglandsale(avglandsale(:,jj)~=0,jj)),...
%         linspace(min(min(avglandsale)),max(max(avglandsale)),10)),10,1);
% end
%     
% h_rent=figure;
% set(h_rent,'Color','white','visible','off');
% for a=1:ERUNS
%     % plot avgrent histogram
%     [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],a);
%     subtightplot(5,5,a,0.01,0.1,0.05);
%     bar(avgrent_mat{ir,ic})
%     axis([0.5 8.5 0 22000]);
%     if isempty(find(ismember([5 10 15 20],a),1))==0
%          set(gca,'xticklabel',[],'yaxislocation','right')
%     elseif a == 25
%         set(gca,'yaxislocation','right')
%     elseif a==23
%         set(gca,'xticklabel',[],'yticklabel',[])
%         xlabel('Avg. Rent t=30 (baseline: (0.1,1.0)','FontSize',14)
%     else
%         set(gca,'xticklabel',[],'yticklabel',[])
%     end
%     if isempty(find(ismember(1:5,a),1))==0
%         title(sprintf('lv(%0.2f)',double(batchparms_unq(a,2))))
%     end
%     icol=ismember([1 6 11 16 21],a);
%     if isempty(find(icol,1))==0
%         ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
%     end 
% end
% 
% h_inc=figure;
% set(h_inc,'visible','off');
% for a=1:ERUNS
%     % plot inc. dist. histogram
%     [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],a);
%     subtightplot(5,5,a,0.01,0.1,0.05);
%     bar(linspace(20000,200000,20),incdist_mat{ir,ic},'histc')
%     axis([10000 210000 0 750])
%     if isempty(find(ismember([5 10 15 20],a),1))==0
%          set(gca,'xticklabel',[],'yaxislocation','right')
%     elseif a == 25
%         set(gca,'yaxislocation','right')
%     elseif a==23
%         set(gca,'xticklabel',[],'yticklabel',[])
%         xlabel('Income t=30 (baseline: (0.1,1.0)','FontSize',14)
%     else
%         set(gca,'xticklabel',[],'yticklabel',[])
%     end
%     if isempty(find(ismember(1:5,a),1))==0
%         title(sprintf('lv(%0.2f)',double(batchparms_unq(a,2))))
%     end
%     icol=ismember([1 6 11 16 21],a);
%     if isempty(find(icol,1))==0
%         ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
%     end 
% end
% saveas(h_inc,'summary_incdist','jpg')
% 
% h_ls=figure;
% set(h_ls,'visible','off');
% for a=1:ERUNS
%     % plot inc. dist. histogram
%     [ic,ir]=ind2sub([EXPTRUNS EXPTRUNS],a);
%     subtightplot(5,5,a,0.01,0.1,0.05);
%     bar(linspace(min(min(avglandsale)),max(max(avglandsale)),10),landsale_mat{ir,ic},'histc')
%     axis([500 27000 0 40])
%     if a == 25
%         set(gca,'yaxislocation','right')
%     elseif a==23
%         set(gca,'xticklabel',[],'yticklabel',[])
%         xlabel('Land Sales t=30 (baseline: (0.1,1.0)','FontSize',14)
%     else
%         set(gca,'xticklabel',[],'yticklabel',[])
%     end
%     if isempty(find(ismember(1:5,a),1))==0
%         title(sprintf('lv(%0.2f)',double(batchparms_unq(a,2))))
%     end
%     icol=ismember([1 6 11 16 21],a);
%     if isempty(find(icol,1))==0
%         ylabel(sprintf('am(%0.3f)',double(batchparms_unq(icol,1))))
%     end 
% end
% saveas(h_ls,'summary_landsales','jpg')
