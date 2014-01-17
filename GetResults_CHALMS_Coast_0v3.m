%%%%%%%% Write results to Excel file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MRUNS=30;
TMAX=30;
TSTART=10;
NWIDTH=80;
NLENGTH=80;
HTYPE=3;
LTYPE=3;
HT=HTYPE*LTYPE;
discount=0.05;

BASELAYER=zeros(NLENGTH,NWIDTH);
COAST=zeros(NLENGTH,NWIDTH);
SCAPE=zeros(NLENGTH,NWIDTH);
COAST(:,1:5)=1;
SCAPE(COAST~=1)=1;
DISTANCE=ones(NLENGTH,NWIDTH);
dist2hznnei=zeros(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
dist2vrtnei=zeros(NLENGTH,NWIDTH);
NZONES=25;
ZONES=zeros(NLENGTH,NWIDTH);
zoning=zeros(NZONES,2);     %[(min lotsize) (max lotsize)]
ZONEMAP=zeros(NLENGTH,NWIDTH);
izoneborder=zeros(1,NWIDTH);
borderdist=zeros(1,NWIDTH);
nzoneslong=sqrt(NZONES);
nzoneswide=sqrt(NZONES);
extralengthz=rem(NLENGTH,nzoneslong);
extrawidthz=rem(NWIDTH,nzoneswide);
izonelength=(length(BASELAYER(:,1))-extralengthz)/nzoneslong;
izonewidth=(length(BASELAYER(1,:))-extrawidthz)/nzoneswide;

zonemarklong=1;
zonemarkwide=1;
for ii=1:nzoneswide
    for jj=1:nzoneslong
        ZONES(zonemarklong:zonemarklong+izonelength-1,...
            zonemarkwide:zonemarkwide+izonewidth-1)=...
            ii*nzoneslong-(nzoneslong-jj);        
        
        if jj==nzoneslong && extralengthz > 0
            ZONES(izonelength*jj+1:izonelength*jj+extralengthz,...
                zonemarkwide:zonemarkwide+izonewidth-1)=...
                ii*nzoneslong-(nzoneslong-jj);        
        end
        if ii==nzoneswide && extrawidthz > 0
            ZONES(zonemarklong:zonemarklong+izonelength-1,...
                izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
                ii*nzoneslong-(nzoneslong-jj);
        end
        zonemarklong=mod(izonelength*jj+1,izonelength*nzoneslong);
    end     
    if jj==nzoneslong && ii==nzoneswide && extralengthz > 0
        ZONES(izonelength*jj+1:izonelength*jj+extralengthz,...
            izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
            ii*nzoneslong-(nzoneslong-jj);
    end
    zonemarkwide=mod(izonewidth*ii+1,izonewidth*nzoneswide);
end

izoned=ismember(ZONES,[3:5 9:10 14:15 19:20 23:25]);
inotzoned=~ismember(ZONES,[3:5 9:10 14:15 19:20 23:25]);
ZONEMAP(izoned)=1;
unzonedacres=length(find(ZONEMAP==0));
zonedacres=length(find(ZONEMAP==1));

northindex=find(ZONEMAP==0);
southindex=find(ZONEMAP==1);

icenterrow=1;
% icentercol=round(NWIDTH/2);
icentercol=6;

for col=1:NWIDTH
    dist2hznnei(1:NLENGTH,col)=abs(col-icentercol).*ones(NLENGTH,1);
end

for row=1:NWIDTH
    dist2vrtnei(row,1:NWIDTH)=abs(row-icenterrow).*ones(1,NLENGTH);
end

for col=1:NWIDTH
    for row=1:NLENGTH
        DISTANCE(row,col)=sqrt(dist2hznnei(row,col)^2+dist2vrtnei(row,col)^2);   
    end
end

for i=1:NWIDTH
    izoneborder(i)=find(ZONEMAP(:,i)==1,1,'first');
    borderdist(i)=DISTANCE(izoneborder(i),i)*0.0395;
end
coastdist=cumsum(SCAPE,2);
coastprox=max(max(coastdist))-coastdist;
%%%%%
MILEDIST=DISTANCE.*0.0395;
%%%%%

testtime=[10 15 20 25 30];
house2acre=[2 0.5 0.2]';
zs=min(max(1./house2acre,0),5);
housesize=[1500 2000 2500]';
zonedist=[0 14 28 42 56 70 84]'.*0.0395;
timeset=(TSTART+1:TMAX);
z=[reshape(repmat(zs,1,3)',HT,1) repmat(housesize,3,1)];

lotsizemap=zeros(NLENGTH,NWIDTH);

fnames=dir;
fnamescell=struct2cell(fnames);
 %%% CHANGE FILE NAME %%%%
h=strmatch('results_CHALMS_Coast_0v3',fnamescell(1,:));

bidsharerlts=zeros(HT,length(TSTART+1:TMAX),length(h));
bidsharestats=zeros(HT,5,length(TSTART+1:TMAX));
bidsharemu=zeros(HT,length(TSTART+1:TMAX));
bidsharesigma=zeros(HT,length(TSTART+1:TMAX));
bidsharemuci=zeros(2,HT,length(TSTART+1:TMAX));
bidsharesigmaci=zeros(2,HT,length(TSTART+1:TMAX));
bidsharese=zeros(HT,length(TSTART+1:TMAX));

htutilvacrlts=zeros(HT,length(TSTART+1:TMAX),length(h));
htutilvacstats=zeros(HT,5,length(TSTART+1:TMAX));
htutilvacmu=zeros(HT,length(TSTART+1:TMAX));
htutilvacsigma=zeros(HT,length(TSTART+1:TMAX));
htutilvacmuci=zeros(2,HT,length(TSTART+1:TMAX));
htutilvacsigmaci=zeros(2,HT,length(TSTART+1:TMAX));
htutilvacse=zeros(HT,length(TSTART+1:TMAX));

htutilnovacrlts=zeros(HT,length(TSTART+1:TMAX),length(h));
htutilnovacstats=zeros(HT,5,length(TSTART+1:TMAX));
htutilnovacmu=zeros(HT,length(TSTART+1:TMAX));
htutilnovacsigma=zeros(HT,length(TSTART+1:TMAX));
htutilnovacmuci=zeros(2,HT,length(TSTART+1:TMAX));
htutilnovacsigmaci=zeros(2,HT,length(TSTART+1:TMAX));
htutilnovacse=zeros(HT,length(TSTART+1:TMAX));

htincomerlts=zeros(HT,length(TSTART+1:TMAX),length(h));
htincomestats=zeros(HT,5,length(TSTART+1:TMAX));
htincomemu=zeros(HT,length(TSTART+1:TMAX));
htincomesigma=zeros(HT,length(TSTART+1:TMAX));
htincomemuci=zeros(2,HT,length(TSTART+1:TMAX));
htincomesigmaci=zeros(2,HT,length(TSTART+1:TMAX));
htincomese=zeros(HT,length(TSTART+1:TMAX));

retrlts=zeros(HT,length(TSTART+1:TMAX),length(h));
retstats=zeros(HT,5,length(TSTART+1:TMAX));
retmu=zeros(HT,length(TSTART+1:TMAX));
retsigma=zeros(HT,length(TSTART+1:TMAX));
retmuci=zeros(2,HT,length(TSTART+1:TMAX));
retsigmaci=zeros(2,HT,length(TSTART+1:TMAX));
retse=zeros(HT,length(TSTART+1:TMAX));

fullrentrlts=zeros(HT,TMAX,length(h));
rentrlts=zeros(HT,length(testtime),length(h));
rentstats=zeros(HT,5,length(testtime));
rentmu=zeros(HT,length(testtime));
rentsigma=zeros(HT,length(testtime));
rentmuci=zeros(2,HT,length(testtime));
rentsigmaci=zeros(2,HT,length(testtime));
rentse=zeros(HT,length(testtime));


zonepctdev=zeros(2,length(testtime));
zonenumlots=zeros(2,length(testtime));
zoneavgprice=zeros(2,length(testtime));
% incomestats=zeros(1,length(testtime));

northrents=zeros(HT,length(testtime),length(h));
nrentstats=zeros(HT,5,length(testtime));
nrentmu=zeros(HT,length(testtime));
nrentsigma=zeros(HT,length(testtime));
nrentmuci=zeros(2,HT,length(testtime));
nrentsigmaci=zeros(2,HT,length(testtime));
nrentse=zeros(HT,length(testtime));

southrents=zeros(HT,length(testtime),length(h));
srentstats=zeros(HT,5,length(testtime));
srentmu=zeros(HT,length(testtime));
srentsigma=zeros(HT,length(testtime));
srentmuci=zeros(2,HT,length(testtime));
srentsigmaci=zeros(2,HT,length(testtime));
srentse=zeros(HT,length(testtime));

vacrlts=zeros(length(testtime),length(h));
pctutilrlts=zeros(length(testtime),length(h));
numlotrlts=zeros(HT,length(testtime),length(h));
lotstats=zeros(HT,5,length(testtime));
lotmu=zeros(HT,length(testtime));
lotsigma=zeros(HT,length(testtime));
lotmuci=zeros(2,HT,length(testtime));
lotsigmaci=zeros(2,HT,length(testtime));
lotse=zeros(HT,length(testtime));

unzonedlotstats=zeros(HT,5,length(testtime));
unzonedlotmu=zeros(HT,length(testtime));
unzonedlotsigma=zeros(HT,length(testtime));
unzonedlotmuci=zeros(2,HT,length(testtime));
unzonedlotsigmaci=zeros(2,HT,length(testtime));
unzonedlotse=zeros(HT,length(testtime));

zonedlotstats=zeros(HT,5,length(testtime));
zonedlotmu=zeros(HT,length(testtime));
zonedlotsigma=zeros(HT,length(testtime));
zonedlotmuci=zeros(2,HT,length(testtime));
zonedlotsigmaci=zeros(2,HT,length(testtime));
zonedlotse=zeros(HT,length(testtime));

lotsumrlts=zeros(length(testtime),length(h));
disprlts=zeros(1,length(h));
maxdistrlts=zeros(1,length(h));
pctdevrlts=zeros(length(testtime),length(h));
testpctdev=zeros(length(testtime),1);
pctbiglotsrlts=zeros(length(testtime),length(h));
plandrlts=zeros(50,4,length(h));
plandinfo=zeros([],5); %[farmerid selltime sellprice dist acres]
avgplandtime=zeros(20,2);
avgplanddist=zeros(50,2);
agretrlts=zeros(2,length(testtime),length(h));
zonedenrlts=zeros(HT,length(h));
zonedenmu=zeros(HT,1);
zonedensigma=zeros(HT,1);
zonedenmuci=zeros(2,HT);
zonedensigmaci=zeros(2,HT);
zonedense=zeros(HT,1);

zoneavgpricestats=zeros(2,5);
overplandstats=zeros(1,5);

unzonedpctdev=zeros(length(testtime),length(h));
zonedpctdev=zeros(length(testtime),length(h));
zonednumlots=zeros(HT,length(testtime),length(h));
unzonednumlots=zeros(HT,length(testtime),length(h));
zoneavgprice=zeros(2,length(h));
overpland=zeros(1,length(h));
incomerlts=zeros(length(testtime),length(h));
outincomerlts=zeros(length(testtime),length(h));

northincome=zeros(length(testtime),length(h));
southincome=zeros(length(testtime),length(h));

% zonedenstats=zeros(HT,length(testtime));
lotsumstats=zeros(length(testtime),5);
lotsummu=zeros(length(testtime),1);
lotsumsigma=zeros(length(testtime),1);
lotsummuci=zeros(2,length(testtime));
lotsumsigmaci=zeros(2,length(testtime));
lotsumse=zeros(length(testtime),1);


vacstats=zeros(length(testtime),5);
vacmu=zeros(length(testtime),1);
vacsigma=zeros(length(testtime),1);
vacmuci=zeros(2,length(testtime));
vacsigmaci=zeros(2,length(testtime));
vacse=zeros(length(testtime),1);

zonedpctdevstats=zeros(length(testtime),5);
zonedpctdevmu=zeros(length(testtime),1);
zonedpctdevsigma=zeros(length(testtime),1);
zonedpctdevmuci=zeros(2,length(testtime));
zonedpctdevsigmaci=zeros(2,length(testtime));
zonedpctdevse=zeros(length(testtime),1);

unzonedpctdevstats=zeros(length(testtime),5);
unzonedpctdevmu=zeros(length(testtime),1);
unzonedpctdevsigma=zeros(length(testtime),1);
unzonedpctdevmuci=zeros(2,length(testtime));
unzonedpctdevsigmaci=zeros(2,length(testtime));
unzonedpctdevse=zeros(length(testtime),1);

incomestats=zeros(length(testtime),5);
incomemu=zeros(length(testtime),1);
incomesigma=zeros(length(testtime),1);
incomemuci=zeros(2,length(testtime));
incomesigmaci=zeros(2,length(testtime));
incomese=zeros(length(testtime),1);

northincomestats=zeros(length(testtime),5);
northincomemu=zeros(length(testtime),1);
northincomesigma=zeros(length(testtime),1);
northincomemuci=zeros(2,length(testtime));
northincomesigmaci=zeros(2,length(testtime));
northincomese=zeros(length(testtime),1);

southincomestats=zeros(length(testtime),5);
southincomemu=zeros(length(testtime),1);
southincomesigma=zeros(length(testtime),1);
southincomemuci=zeros(2,length(testtime));
southincomesigmaci=zeros(2,length(testtime));
southincomese=zeros(length(testtime),1);

outincomestats=zeros(length(testtime),5);
outincomemu=zeros(length(testtime),1);
outincomesigma=zeros(length(testtime),1);
outincomemuci=zeros(2,length(testtime));
outincomesigmaci=zeros(2,length(testtime));
outincomese=zeros(length(testtime),1);

utilstats=zeros(length(testtime),5);
utilmu=zeros(length(testtime),1);
utilsigma=zeros(length(testtime),1);
utilmuci=zeros(2,length(testtime));
utilsigmaci=zeros(2,length(testtime));
utilse=zeros(length(testtime),1);

devstats=zeros(length(testtime),5);
devmu=zeros(length(testtime),1);
devsigma=zeros(length(testtime),1);
devmuci=zeros(2,length(testtime));
devsigmaci=zeros(2,length(testtime));
devse=zeros(length(testtime),1);

biglotsstats=zeros(length(testtime),5);
biglotsmu=zeros(length(testtime),1);
biglotssigma=zeros(length(testtime),1);
biglotsmuci=zeros(2,length(testtime));
biglotssigmaci=zeros(2,length(testtime));
biglotsse=zeros(length(testtime),1);

agretsoldstats=zeros(length(testtime),5);
agretsoldmu=zeros(length(testtime),1);
agretsoldsigma=zeros(length(testtime),1);
agretsoldmuci=zeros(2,length(testtime));
agretsoldsigmaci=zeros(2,length(testtime));
agretsoldse=zeros(length(testtime),1);

agretfarmstats=zeros(length(testtime),5);
agretfarmmu=zeros(length(testtime),1);
agretfarmsigma=zeros(length(testtime),1);
agretfarmmuci=zeros(2,length(testtime));
agretfarmsigmaci=zeros(2,length(testtime));
agretfarmse=zeros(length(testtime),1);

Exptretrlts=zeros(HT,TMAX,length(h));
Realretrlts=zeros(HT,TMAX,length(h));
Nlots=zeros(HT,TMAX,length(h));
IDEALSET=zeros(HT,length(testtime),length(h));
PROFSET=zeros(HT,length(testtime),length(h));
AVGBVAR=zeros(HT,TMAX,length(h));

lottypemap=zeros(NLENGTH,NWIDTH,HT,length(testtime));
secprobmap=zeros(NLENGTH,NWIDTH,HT,length(testtime));
subvarmap=zeros(NLENGTH,NWIDTH,length(testtime));

VARMAP=zeros(NLENGTH,NWIDTH,length(testtime));

VARLAYERSAVE=zeros(NLENGTH,NWIDTH);
avgthresh=[1 0.36 0.4 0.4 0.45];
% avgthresh=[1 0.33 0.4 0.47 0.63];
AVGMAP=zeros(NLENGTH,NWIDTH,length(testtime));
CONFMAP=zeros(NLENGTH,NWIDTH,length(testtime));
carcost=zeros(TMAX,length(h));
carrycostrlts=zeros(1,TMAX);

PMAPS=zeros(NLENGTH,NWIDTH,length(testtime),length(h));
avgPREFMAP=zeros(NLENGTH,NWIDTH,length(testtime));

t30maps=zeros(NLENGTH,NWIDTH,length(h));

%%%%%%
avglotsize=zeros(35,length(h));
plandint=zeros(size(avglotsize));
pctdevdist=zeros(size(avglotsize));
houseperacre=zeros(size(avglotsize));
meandispersion=zeros(length(h),1);
plandzone=zeros(2,length(timeset));

fullnumlots=zeros(HT,length(TSTART:TMAX),length(h));
%%%%%%

for mr=1:length(h)
    filename=char(fnamescell(1,h(mr)));
    load(eval('filename'))
    
    %%%%%% Density and Pct Dev w/ distance %%%%%%%%%%%%%%%%%5
    subRlottype=Rlottype(:,:,length(testtime));
    subdensity=subRlottype;
    ifinddev=find(subRlottype ~= 0);
    subdensity(ifinddev)=1./z(subRlottype(ifinddev),1);
    
    startpt=0;
%     startpt=0.5195;
    endpt=startpt+0.1;
    
%     while startpt < max(max(MILEDIST))
    for ii=1:35
        iring=find(MILEDIST >= startpt & MILEDIST < endpt);
        idev=(subRlottype(iring)~=0);
        intfarm=unique(setupmap(iring));
        ifarm=ismember(plandinfo(:,1),intfarm);
        plandint(ii,mr)=sum(plandinfo(ifarm,3).*plandinfo(ifarm,5))/sum(...
            plandinfo(ifarm,5));
        pctdevdist(ii,1,mr)=length(find(idev==1))/length(iring);
        houseperacre(ii,mr)=sum(1./z(subRlottype(iring(idev))))/length(iring);
%         for lt=1:HT
%         lotden(lt)=length(find(subRlottype(iring(idev))==lt))/z(lt,1);
%         end
        avglotsize(ii,mr)=1/(sum(1./z(subRlottype(iring(idev))))/length(iring(idev)));
%         avglotsize(ii,1,mr)=mean(z(subRlottype(iring(idev)),1));
%         avglotden(ii,1,mr)=mean(lotden);
        startpt=startpt+0.1;
        endpt=startpt+0.1;
        if startpt >= max(max(MILEDIST(ifinddev)))
            break
        end
    end
    finedistances=(0.1:0.1:endpt-0.1);
%     subavglotsize=zeros(size(avglotsize));
%     subavglotsize(~isnan(avglotsize))=avglotsize(~isnan(avglotsize));
%     meanlotsize=mean(subavglotsize,3)';
%     
%     subhouseperacre=zeros(size(houseperacre));
%     subhouseperacre(~isnan(houseperacre))=houseperacre(~isnan(houseperacre));
%     meanhouseperacre=mean(subhouseperacre,3)';
    %%%%%% Mean Dispersion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALTDISTANCE=ones(NLENGTH,NWIDTH);
    altdist2hznnei=zeros(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
    altdist2vrtnei=zeros(NLENGTH,NWIDTH);

    %%% Conversion to kilometers to match Irwin and Bockstael (2007)
    %%% neighborhood of 1 km^2
    %%% 1kilometer = 0.6214 mi
    %%% Radius of 1 km^2 area = sqrt(2)/2 = 0.7071
    %%% Dispersion length = 0.7071*0.6214 = 0.4394
    neilength=0.4394;

    dispersion=zeros(length(find(ifinddev)),1);
    for ic=1:length(ifinddev)

        [icenterrow,icentercol]=ind2sub([NLENGTH NWIDTH],ifinddev(ic));
        altdist2hznnei=ones(NLENGTH,1)*abs((1:NWIDTH)-icentercol);
        altdist2vrtnei=abs((1:NLENGTH)-icenterrow)'*ones(1,NLENGTH);
        ALTDISTANCE=sqrt(altdist2hznnei.^2+altdist2vrtnei.^2);
%         for col=1:NWIDTH
%         for col=max(1,icentercol-15):min(NWIDTH,icentercol+15)
%             altdist2hznnei(1:NLENGTH,col)=abs(col-icentercol).*ones(NLENGTH,1);
%         end
% 
% %         for row=1:NWIDTH
%         for row=max(1,icentercol-15):min(NWIDTH,icentercol+15)
%             altdist2vrtnei(row,1:NWIDTH)=abs(row-icenterrow).*ones(1,NLENGTH);
%         end

%         for col=1:NWIDTH
%         for col=max(1,icentercol-15):min(NWIDTH,icentercol+15)
% %             for row=1:NLENGTH
%             for row=max(1,icentercol-15):min(NWIDTH,icentercol+15)
%                 ALTDISTANCE(row,col)=sqrt(altdist2hznnei(row,col)^2+altdist2vrtnei(row,col)^2);
%             end
%         end
        ALTDISTANCE=ALTDISTANCE.*0.0395; % convert to miles
        
        ineihood=(ALTDISTANCE <= neilength & ALTDISTANCE > 0);

        dispersion(ic)=length(find(subRlottype(ineihood)==0))/...
            length(find(ineihood==1));
    end
    
    meandispersion(mr)=sum(dispersion)/length(ifinddev);
    
    for itt=TSTART:TMAX
        ifullrents=~isnan(avgrent(:,itt));
        fullrentrlts(ifullrents,itt,mr)=avgrent(ifullrents,itt);
    end
    
    for itt2=1:20
        iretrlts=~isnan(Realavgret(:,TSTART+itt2));
        retrlts(iretrlts,itt2,mr)=Realavgret(iretrlts,TSTART+itt2);
        
%         ihtutilvacrlts=~isnan(htutilwvac(:,TSTART+itt2));
%         htutilvacrlts(ihtutilvacrlts,itt2,mr)=htutilwvac(ihtutilvacrlts,TSTART+itt2);
%         
%         ihtutilnovacrlts=~isnan(htutilnovac(:,TSTART+itt2));
%         htutilnovacrlts(ihtutilnovacrlts,itt2,mr)=htutilnovac(ihtutilnovacrlts,TSTART+itt2);
    
        ihtincomerlts=~isnan(htincome(:,TSTART+itt2));
        htincomerlts(ihtincomerlts,itt2,mr)=htincome(ihtincomerlts,TSTART+itt2);
    
        ibidsharerlts=~isnan(bidshare(:,TSTART+itt2));
        bidsharerlts(ibidsharerlts,itt2,mr)=bidshare(ibidsharerlts,TSTART+itt2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    inzoned=unique(setupmap(izoned));
    outzoned=unique(setupmap(inotzoned));
    outzoned=outzoned(outzoned~=0);
    
    t30maps(:,:,mr)=Rlottype(:,:,length(testtime));
    
    for tt=1:length(testtime)
        subzonelots=Rlottype(:,:,tt);
        subzonerents=Rrent(:,:,tt);
        subincome=Rincome(:,:,tt);
    
        irents=~isnan(avgrent(:,testtime(tt)));
        rentrlts(irents,tt,mr)=avgrent(irents,testtime(tt));
        vacrlts(tt,mr)=mean(vacrate(TSTART:testtime(tt)));
%         pctutilrlts(tt,mr)=mean(pctutildiff(TSTART:testtime(tt)));
        pctdevrlts(tt,mr)=length(find(Rbaselayer(:,:,tt)==1))/(NLENGTH*NWIDTH);
        pctbiglotsrlts(tt,mr)=length(find(Rlottype(:,:,tt) > 9))/...
            length(find(Rbaselayer(:,:,tt)==1));

        subbaselayer=Rbaselayer(:,:,tt);
        unzonedpctdev(tt,mr)=length(find(subbaselayer(inotzoned)==1))/unzonedacres;
        zonedpctdev(tt,mr)=length(find(subbaselayer(izoned)==1))/zonedacres;
        incomerlts(tt,mr)=consumerstats(2,testtime(tt));
        outincomerlts(tt,mr)=oldincome(testtime(tt));
%         farmerset=(1:50)';
%         idfarmnotsold=~ismember(farmerset,farmsoldinfo(:,1));
%         ifarmnotsold=ismember(setupmap,farmerset(idfarmnotsold));
        ifarmnotsold=unique(setupmap(Rbaselayer(:,:,tt)==0));
        ifarmnotsold=ifarmnotsold(ifarmnotsold~=0);
        ifarmsold=unique(setupmap(Rbaselayer(:,:,tt)==1));
        ifarmsold=ifarmsold(ifarmsold~=0);
        subpland=Rpland(:,:,1);
        agretrlts(1,tt,mr)=mean(subpland(ismember(setupmap,ifarmsold)));
        agretrlts(2,tt,mr)=mean(subpland(ismember(setupmap,ifarmnotsold)));
%         agretrlts(1,tt,mr)=mean(Farmstats(ifarmsold,1,1));
%         agretrlts(2,tt,mr)=mean(Farmstats(ifarmnotsold,1,1));
        subcountmap=subvarmap(:,:,tt);
        ivarlayer=(Rbaselayer(:,:,tt)==1);
        subcountmap(ivarlayer)=subcountmap(ivarlayer)+1;
        subvarmap(:,:,tt)=subcountmap;

        for lt=1:HT
            ilots=find(subzonelots==lt);
            inorthlots=ismember(ilots,northindex);
            if isempty(find(inorthlots,1))==0
                if tt==1
                    northrents(lt,tt,mr)=avgrent(lt,testtime(tt));
                else
                    northrents(lt,tt,mr)=mean(subzonerents(ilots(inorthlots)));
                end
            end
            isouthlots=ismember(ilots,southindex);
            if isempty(find(isouthlots,1))==0
                southrents(lt,tt,mr)=mean(subzonerents(ilots(isouthlots)));
            end
            
            iltsize=find(subzonelots == lt);                                     
            lotsizemap(iltsize)=z(lt,1);
            
            numlotrlts(lt,tt,mr)=length(find(Rlottype(:,:,tt)==lt))/z(lt,1);
            sublotsmap=Rlottype(:,:,tt);
            unzonednumlots(lt,tt,mr)=length(find(sublotsmap(inotzoned)==lt))/z(lt,1);
            zonednumlots(lt,tt,mr)=length(find(sublotsmap(izoned)==lt))/z(lt,1);
            sublottypemap=lottypemap(:,:,lt,tt);
            ilottypemap=(Rlottype(:,:,tt)==lt);
            sublottypemap(ilottypemap)=sublottypemap(ilottypemap)+1;
            lottypemap(:,:,lt,tt)=sublottypemap;  
        end
        lotsumrlts(tt,mr)=sum(numlotrlts(:,tt,mr));
        fullnumlots(:,1,mr)=numlotrlts(:,1,mr);
%         irents=~isnan(avgrent(:,testtime(tt)));
%         rentrlts(irents,tt,mr)=sum(avgrent(irents,testtime(tt)).*...
%             numlotrlts(irents,tt,mr))./lotsumrlts(tt,mr);
        if tt >=2
            inorthincome=subincome(inotzoned);
            inorthltsize=lotsizemap(inotzoned);
            northincome(tt,mr)=mean(sum(inorthincome(inorthincome~=0).*...
                (1./inorthltsize(inorthincome~=0)))/sum(1./inorthltsize(inorthincome~=0)));
            isouthincome=subincome(izoned);
            isouthltsize=lotsizemap(izoned);
            southincome(tt,mr)=mean(sum(isouthincome(isouthincome~=0).*...
                (1./isouthltsize(isouthincome~=0)))/sum(1./isouthltsize(isouthincome~=0)));
        end
        
        PMAPS(:,:,tt,mr)=PREFMAP(:,:,tt);
        
    end
    totacres=length(find(ismember(setupmap,farmsoldinfo(:,1))==1));
    
    zoneavgprice(1,mr)=mean(farmsoldinfo(ismember(farmsoldinfo(:,1),outzoned),4));
    zoneavgprice(2,mr)=mean(farmsoldinfo(ismember(farmsoldinfo(:,1),inzoned),4));

%     acrewght=zeros([],1);
    nfacres=zeros([],1);
    for nf=1:length(farmsoldinfo(:,1))
       nfacres(nf,1)=length(find(ismember(setupmap,farmsoldinfo(nf,1))==1)); 
    end
    plandrlts(farmsoldinfo(:,1),:,mr)=farmsoldinfo(:,[1 2 4 7]);
    plandinfo(length(plandinfo(:,1))+1:length(plandinfo(:,1))+...
        length(farmsoldinfo(:,1)),:)=...
        [farmsoldinfo(:,[1 2 4 7]) nfacres];
    
    %Overall, weigthed average land price
    overpland(mr)=sum(plandinfo(:,3).*plandinfo(:,5))/sum(plandinfo(:,5));
    
    fullnumlots(:,2:length(TSTART:TMAX),mr)=fullnumlots(:,1,mr)*...
        ones(1,length(TSTART+1:TMAX))+cumsum(htperyear(:,TSTART+1:TMAX),2);
    
    disprlts(mr)=meandispersion(mr);
    maxdistrlts(mr)=maxdevdist;
    zonedenrlts(:,mr)=zonedensity;
    
    Exptretrlts(:,:,mr)=Exptret/discount;
    Realretrlts(:,:,mr)=Realavgret;
    IDEALSET(:,:,mr)=idealset;
    PROFSET(:,:,mr)=profset;
    AVGBVAR(:,:,mr)=avgbrokervar;
    carcost(:,mr)=carrycost;
    
    clear('consumerstats','vacstats','BT','Rvacland','Rrent','Rreturn',...
        'Rlottype','Rincome','Rbaselayer','Rpop','Rpland','Rvacrate','Rvaclots',...
        'Rnumlot','Rleftoverpop','avgrentdynms','rentdynstats','deltalots',...
        'farmsoldinfo','avgrent','avgfarmsize','stdfarmsize','DELTA',...
        'survivalrate','LOCWGHT','REGWGHT','PCTSEARCH','zeta','HIBETA',...
        'MIDBETA','LOWBETA','POPGROW','ccost','newhouseset','bidtot',...
        'meandisp','maxdevdist','utilgini','incgini','zonedensity',...
        'pctutildiff','VARLAYER','vacrate','Farmstats','oldincome','Realreturn',...
        'Realavgret','Exptrentdiff','Avgexptdiff','Newlottype','htincome',...
        'numtotbids','totfarmrecord','htperyear','Newbidlevel','Dynltmap',...
        'Dynrentmap','totbrokerrecord','Farmdist2dev','Bprojmap','Bmodelmap',...
        'Dynplandmap','WTAlandmap','WTPlandmap','Landprojmap','Landmodmap',...
        'bidshare','buildshare','farmacreinfo','totfarmrecord','landdemand',...
        'Exptprofit','Exptret','Realexptret','Realavgexptret','idealset',...
        'profset','avgbrokervar','carrycost','con2lot','CONINFO','PREFMAP')
end

%%%%%% Build VARMAP for each time mark %%%%%
agretrlts(1,1,:)=0;
for tt=1:length(testtime)
    VARMAP(:,:,tt)=subvarmap(:,:,tt)./length(h);



    %%%%%% Load last run's file %%%%%%%%%
    % load resultsLOC30_160210_20.mat

    % VARLAYERSAVE(:,:)=VARLAYER;
    % varmap=VARLAYERSAVE/max(max(VARLAYER));
    ivarthresh=(VARMAP(:,:,tt) >= avgthresh(tt));
    locmat=find(ivarthresh==1);

    lottypemap(:,:,:,tt)=lottypemap(:,:,:,tt)./length(h);
    [maxprob,imaxprob]=max(lottypemap(:,:,:,tt),[],3);


    for lt=1:HT
        singleltmap=lottypemap(:,:,lt,tt);
        iltzero=(imaxprob==lt);
        singleltmap(iltzero)=0;
        secprobmap(:,:,lt,tt)=singleltmap;
    end
    [secprob,isecprob]=max(secprobmap(:,:,:,tt),[],3);

    confprob=maxprob-secprob;
%     PRODMAP=Rpland(:,:,1);
    %%%%% Test representativeness of varmap
    testpctdev(tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
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
    if tt > 1
        idevthere=(AVGMAP(:,:,tt-1)~=0);
        indthresh=find(ivarthresh==1);
        inddev=find(idevthere==1);
        idevadd=~ismember(indthresh,inddev);
        subAVGMAP=AVGMAP(:,:,tt-1);
        subAVGMAP(indthresh(idevadd))=imaxprob(indthresh(idevadd));
        AVGMAP(:,:,tt)=subAVGMAP;
    else
        subAVGMAP(ivarthresh)=imaxprob(ivarthresh);
        AVGMAP(:,:,tt)=subAVGMAP;
    end
    testmeandisp=sum(altprop)/length(imapcover);
    subCONFMAP=zeros(NLENGTH,NWIDTH);
    subCONFMAP(ivarthresh)=confprob(ivarthresh);
    CONFMAP(:,:,tt)=subCONFMAP;
    
    timePMAPS=zeros(NLENGTH,NWIDTH,length(h));
    timePMAPS(:,:,:)=PMAPS(:,:,tt,:);
    for il=1:length(locmat)
        [irow,icol]=ind2sub([NLENGTH NWIDTH],locmat(il));
        isnotzero=(timePMAPS(irow,icol,:)~=0);
        avgPREFMAP(irow,icol,tt)=mean(timePMAPS(irow,icol,isnotzero));
    end
    ifindnan=find(isnan(avgPREFMAP(:,:,tt))==1);  
end

Exptretrlts(isnan(Exptretrlts))=0;
Realretrlts(isnan(Realretrlts))=0;
DevProf=numlotrlts.*Realretrlts(:,testtime,:);
ER=zeros(HT,TMAX-TSTART);
RR=zeros(HT,TMAX-TSTART);

for it=TSTART+1:TMAX
    for lt=1:HT
        subER=Exptretrlts(lt,it,:);
        ER(lt,it-TSTART)=mean(subER(subER~=0));
        subRR=Realretrlts(lt,it,:);
        RR(lt,it-TSTART)=mean(subRR(subRR~=0));
    end
end
ER(isnan(ER))=0;
RR(isnan(RR))=0;
Avgexptret=mean(ER,2);
Avgrealret=mean(RR,2);
newlots=numlotrlts(:,1,:);
newlots(:,2:5,:)=diff(numlotrlts,1,2);
% totdevprof=sum(RR,2);
totdevprof=sum(mean(newlots,3).*RR(:,[1 5 10 15 20]),2);
totexptprof=sum(mean(newlots,3).*ER(:,[1 5 10 15 20]),2);

% offerset(:,1:30)=numlotrlts(:,5,:);
offerset=mean(numlotrlts,3);
avgidealset=mean(IDEALSET,3);
diffset=mean(IDEALSET-numlotrlts,3);

avgprofset=mean(PROFSET,3);
diffprofset=mean(PROFSET-numlotrlts,3);

avgvar=mean(AVGBVAR(:,10:30,:),3);

%%
    
figure(1)
surf(AVGMAP(:,:,tt))
axis ij
view(0,90)
title('Average Landscape Outcome, t=30')
set(gca,'clim',[0 HT])
colorbar

figure(2)
surf(VARMAP(:,:,5))
axis ij
view(0,90)
title('Probability of Development, t=30')
colorbar

figure(3)
pcolor(avgPREFMAP(:,:,tt))
axis ij
title('Average Consumer Preference for Proximity to Coast')
set(gca,'clim',[0 0.3])
colorbar

% % 
% % figure(2)
% surf(CONFMAP(:,:,2))
% axis ij
% view(0,90)
% title('Likelihood of Average Landscape Outcome, No Zoning, t=15')
% colorbar


%% Generate Statistics %%

subavglotsize=zeros(size(avglotsize));
subavglotsize(~isnan(avglotsize))=avglotsize(~isnan(avglotsize));
subhouseperacre=zeros(size(houseperacre));
subhouseperacre(~isnan(houseperacre))=houseperacre(~isnan(houseperacre));
subplandint=zeros(size(avglotsize));
subplandint(~isnan(plandint))=plandint(~isnan(plandint));

northfarmers=unique(setupmap(ZONEMAP==0));
northfarmers=northfarmers(northfarmers~=0);
southfarmers=unique(setupmap(ZONEMAP==1));
inorthfarmers=ismember(plandinfo(:,1),northfarmers);
isouthfarmers=ismember(plandinfo(:,1),southfarmers);


for imr=1:35
    inonzerolot=(subavglotsize(imr,:) ~= 0);
    inonzerohouse=(subhouseperacre(imr,:) ~= 0);
    inonzeroplandint=(subplandint(imr,:) ~= 0);
    meanlotsize(imr)=mean(subavglotsize(imr,inonzerolot));
    meanhouseperacre(imr)=mean(subhouseperacre(imr,inonzerohouse));
    meanplandint(imr)=mean(subplandint(imr,inonzeroplandint));
end

planddist=unique(plandinfo(:,4));


for pt=1:20
   avgplandtime(pt,1)=timeset(pt);
   %%% weight by acreage
   itimeset=(plandinfo(:,2)==timeset(pt));
   infarm=inorthfarmers.*itimeset;
   isfarm=isouthfarmers.*itimeset;
   plandzone(1,pt)=mean(plandinfo(infarm==1,3));
   plandzone(2,pt)=mean(plandinfo(isfarm==1,3));
   
%    avgplandtime(pt,2)=mean(plandinfo(itimeset,3).*plandinfo(itimeset,5));
   avgplandtime(pt,2)=sum(plandinfo(itimeset,3).*plandinfo(itimeset,5))/...
       sum(plandinfo(itimeset,5));
   
   for lt=1:HT
       subrets(1:length(h))=retrlts(lt,pt,:);
       iposrets=(subrets~=0);
%        subhtvac(1:length(h))=htutilvacrlts(lt,pt,:);
%        iposhtv=(subhtvac~=0);
%        subhtnovac(1:length(h))=htutilnovacrlts(lt,pt,:);
%        iposhtnv=(subhtnovac~=0);
       subhtinc(1:length(h))=htincomerlts(lt,pt,:);
       iposhti=(subhtinc~=0);
       subbids(1:length(h))=bidsharerlts(lt,pt,:);
       iposbids=(subbids~=0);
       
       subretlots(1:length(h))=fullnumlots(lt,pt+1,:);
       iposretlots=(subretlots~=0);
       
       if isempty(find(iposrets,1))==1
           continue
       end
       
       [retmu(lt,pt),retsigma(lt,pt),retmuci(:,lt,pt),retsigmaci(:,lt,pt)]=normfit(subrets(iposrets));
       retmu(lt,pt)=sum(subrets(iposrets).*subretlots(iposrets))/sum(subretlots(iposrets));
       retse(lt,pt)=retsigma(lt,pt)/sqrt(length(subrets(iposrets)));
       
%        [htutilvacmu(lt,pt),htutilvacsigma(lt,pt),htutilvacmuci(:,lt,pt),...
%            htutilvacsigmaci(:,lt,pt)]=normfit(subhtvac(iposhtv));
%        htutilvacmu(lt,pt)=sum(subhtvac(iposhtv).*subretlots(iposhtv))/sum(subretlots(iposhtv));
%        htutilvacse(lt,pt)=htutilvacsigma(lt,pt)/sqrt(length(subhtvac(iposhtv)));
%        
%        [htutilnovacmu(lt,pt),htutilnovacsigma(lt,pt),htutilnovacmuci(:,lt,pt),...
%            htutilnovacsigmaci(:,lt,pt)]=normfit(subhtnovac(iposhtnv));
%        htutilnovacmu(lt,pt)=sum(subhtnovac(iposhtnv).*subretlots(iposhtnv))/sum(subretlots(iposhtnv));
%        htutilnovacse(lt,pt)=htutilnovacsigma(lt,pt)/sqrt(length(subhtnovac(iposhtnv)));
       
       [htincomemu(lt,pt),htincomesigma(lt,pt),htincomemuci(:,lt,pt),...
           htincomesigmaci(:,lt,pt)]=normfit(subhtinc(iposhti));
       htincomemu(lt,pt)=sum(subhtinc(iposhti).*subretlots(iposhti))/sum(subretlots(iposhti));
       htincomese(lt,pt)=htincomesigma(lt,pt)/sqrt(length(subhtinc(iposhti)));
       
       [bidsharemu(lt,pt),bidsharesigma(lt,pt),bidsharemuci(:,lt,pt),...
           bidsharesigmaci(:,lt,pt)]=normfit(subbids(iposbids));
       bidsharemu(lt,pt)=sum(subbids(iposbids).*subretlots(iposbids))/sum(subretlots(iposbids));
       bidsharese(lt,pt)=bidsharesigma(lt,pt)/sqrt(length(subbids(iposbids)));
   end
   retstats(:,:,pt)=max([retmu(:,pt),retsigma(:,pt),retmuci(1,:,pt)',retmuci(2,:,pt)',retse(:,pt)],0);
   htutilvacstats(:,:,pt)=max([htutilvacmu(:,pt),htutilvacsigma(:,pt),...
       htutilvacmuci(1,:,pt)',htutilvacmuci(2,:,pt)',htutilvacse(:,pt)],0);
   htutilnovacstats(:,:,pt)=max([htutilnovacmu(:,pt),htutilnovacsigma(:,pt),...
       htutilnovacmuci(1,:,pt)',htutilnovacmuci(2,:,pt)',htutilnovacse(:,pt)],0);
   htincomestats(:,:,pt)=max([htincomemu(:,pt),htincomesigma(:,pt),...
       htincomemuci(1,:,pt)',htincomemuci(2,:,pt)',htincomese(:,pt)],0);
   bidsharestats(:,:,pt)=max([bidsharemu(:,pt),bidsharesigma(:,pt),...
       bidsharemuci(1,:,pt)',bidsharemuci(2,:,pt)',bidsharese(:,pt)],0);
end
plandzone(isnan(plandzone))=0;
% for pd=1:50
for pd=1:length(planddist)
   avgplanddist(pd,1)=planddist(pd).*0.0395;
   avgplanddist(pd,2)=mean(plandinfo(plandinfo(:,4)==planddist(pd),3));
end

carrycostrlts(1,:)=mean(carcost,2)';
for tt=1:length(testtime)
    for lt=1:HT
        subnorthrents(1:length(h))=northrents(lt,tt,:);
        iposnorthrents=(subnorthrents~=0);
        subsouthrents(1:length(h))=southrents(lt,tt,:);
        ipossouthrents=(subsouthrents~=0);
        subnorthlots(1:length(h))=unzonednumlots(lt,tt,:);
        iposnorthlots=(subnorthlots~=0);
        subsouthlots(1:length(h))=zonednumlots(lt,tt,:);
        ipossouthlots=(subsouthlots~=0);
        subrents(1:length(h))=rentrlts(lt,tt,:);
        iposrents=(subrents~=0);
        sublots(1:length(h))=numlotrlts(lt,tt,:);
        iposlots=(sublots~=0);
        
        [nrentmu(lt,tt),nrentsigma(lt,tt),nrentmuci(:,lt,tt),nrentsigmaci(:,lt,tt)]=normfit(subnorthrents(iposnorthrents));
        nrentmu(lt,tt)=sum(subnorthrents(iposnorthrents).*subnorthlots(iposnorthrents))/sum(subnorthlots(iposnorthrents));
        nrentse(lt,tt)=nrentsigma(lt,tt)/sqrt(length(subnorthrents(iposnorthrents)));
        
        [srentmu(lt,tt),srentsigma(lt,tt),srentmuci(:,lt,tt),srentsigmaci(:,lt,tt)]=normfit(subsouthrents(ipossouthrents));
        srentmu(lt,tt)=sum(subsouthrents(ipossouthrents).*subsouthlots(ipossouthrents))/sum(subsouthlots(ipossouthrents));
        srentse(lt,tt)=srentsigma(lt,tt)/sqrt(length(subsouthrents(ipossouthrents)));
        
        if isempty(find(iposrents,1))==1
            continue
        end
        [rentmu(lt,tt),rentsigma(lt,tt),rentmuci(:,lt,tt),rentsigmaci(:,lt,tt)]=normfit(subrents(iposrents));
        rentmu(lt,tt)=sum(subrents(iposrents).*sublots(iposrents))/sum(sublots(iposrents));
        rentse(lt,tt)=rentsigma(lt,tt)/sqrt(length(subrents(iposrents)));
%         [lotmu(lt,tt),lotsigma(lt,tt),lotmuci(:,lt,tt),lotsigmaci(:,lt,tt)]=normfit(sublots(iposlots));
        [lotmu(lt,tt),lotsigma(lt,tt),lotmuci(:,lt,tt),lotsigmaci(:,lt,tt)]=normfit(sublots);        
        lotse(lt,tt)=lotsigma(lt,tt)/sqrt(length(sublots));
        subunzonednumlots(1:length(h))=unzonednumlots(lt,tt,:);
        subzonednumlots(1:length(h))=zonednumlots(lt,tt,:);
        [unzonedlotmu(lt,tt),unzonedlotsigma(lt,tt),unzonedlotmuci(:,lt,tt),unzonedlotsigmaci(:,lt,tt)]=normfit(subunzonednumlots);        
        unzonedlotse(lt,tt)=unzonedlotsigma(lt,tt)/sqrt(length(subunzonednumlots));
        [zonedlotmu(lt,tt),zonedlotsigma(lt,tt),zonedlotmuci(:,lt,tt),zonedlotsigmaci(:,lt,tt)]=normfit(subzonednumlots);        
        zonedlotse(lt,tt)=zonedlotsigma(lt,tt)/sqrt(length(subzonednumlots));
        
    end
    for zl=1:length(zonedenrlts(:,1))
        iposzone=~isnan(zonedenrlts(zl,:));
        [zonedenmu(zl),zonedensigma(zl),zonedenmuci(:,zl),zonedensigmaci(:,zl)]=normfit(zonedenrlts(zl,iposzone));
        zonedense(zl)=zonedensigma(zl)/sqrt(length(h));
    end
    

    % keyboard    %may need to transpose, check size of rentmu
    nrentstats(:,:,tt)=max([nrentmu(:,tt),nrentsigma(:,tt),nrentmuci(1,:,tt)',nrentmuci(2,:,tt)',nrentse(:,tt)],0);
    srentstats(:,:,tt)=max([srentmu(:,tt),srentsigma(:,tt),srentmuci(1,:,tt)',srentmuci(2,:,tt)',srentse(:,tt)],0);
        
    rentstats(:,:,tt)=max([rentmu(:,tt),rentsigma(:,tt),rentmuci(1,:,tt)',rentmuci(2,:,tt)',rentse(:,tt)],0);
    lotstats(:,:,tt)=max([lotmu(:,tt),lotsigma(:,tt),lotmuci(1,:,tt)',lotmuci(2,:,tt)',lotse(:,tt)],0);
    
    
    unzonedlotstats(:,:,tt)=max([unzonedlotmu(:,tt),unzonedlotsigma(:,tt),unzonedlotmuci(1,:,tt)',unzonedlotmuci(2,:,tt)',unzonedlotse(:,tt)],0);
    zonedlotstats(:,:,tt)=max([zonedlotmu(:,tt),zonedlotsigma(:,tt),zonedlotmuci(1,:,tt)',zonedlotmuci(2,:,tt)',zonedlotse(:,tt)],0);
    
    [lotsummu(tt),lotsumsigma(tt),lotsummuci(:,tt),lotsumsigmaci(:,tt)]=normfit(lotsumrlts(tt,:));
    lotsumse(tt)=lotsumsigma(tt)/sqrt(length(h));
    lotsumstats(tt,:)=max([lotsummu(tt),lotsumsigma(tt),lotsummuci(1,tt),lotsummuci(2,tt),lotsumse(tt)],0);

    [vacmu(tt),vacsigma(tt),vacmuci(:,tt),vacsigmaci(:,tt)]=normfit(vacrlts(tt,:));
    vacse(tt)=vacsigma(tt)/sqrt(length(h));
    vacstats(tt,:)=max([vacmu(tt),vacsigma(tt),vacmuci(1,tt),vacmuci(2,tt),vacse(tt)],0);
    
    [incomemu(tt),incomesigma(tt),incomemuci(:,tt),incomesigmaci(:,tt)]=normfit(incomerlts(tt,:));
    incomese(tt)=incomesigma(tt)/sqrt(length(h));
    incomestats(tt,:)=max([incomemu(tt),incomesigma(tt),incomemuci(1,tt),incomemuci(2,tt),incomese(tt)],0);
    
    [northincomemu(tt),northincomesigma(tt),northincomemuci(:,tt),northincomesigmaci(:,tt)]=normfit(northincome(tt,:));
    northincomese(tt)=northincomesigma(tt)/sqrt(length(h));
    northincomestats(tt,:)=max([northincomemu(tt),northincomesigma(tt),northincomemuci(1,tt),northincomemuci(2,tt),northincomese(tt)],0);
    
    [southincomemu(tt),southincomesigma(tt),southincomemuci(:,tt),southincomesigmaci(:,tt)]=normfit(southincome(tt,:));
    southincomese(tt)=southincomesigma(tt)/sqrt(length(h));
    southincomestats(tt,:)=max([southincomemu(tt),southincomesigma(tt),southincomemuci(1,tt),southincomemuci(2,tt),southincomese(tt)],0);
    
    [outincomemu(tt),outincomesigma(tt),outincomemuci(:,tt),outincomesigmaci(:,tt)]=normfit(outincomerlts(tt,:));
    outincomese(tt)=outincomesigma(tt)/sqrt(length(h));
    outincomestats(tt,:)=max([outincomemu(tt),outincomesigma(tt),outincomemuci(1,tt),outincomemuci(2,tt),outincomese(tt)],0);
    
    [zonedpctdevmu(tt),zonedpctdevsigma(tt),zonedpctdevmuci(:,tt),zonedpctdevsigmaci(:,tt)]=normfit(zonedpctdev(tt,:));
    zonedpctdevse(tt)=zonedpctdevsigma(tt)/sqrt(length(h));
    zonedpctdevstats(tt,:)=max([zonedpctdevmu(tt),zonedpctdevsigma(tt),zonedpctdevmuci(1,tt),zonedpctdevmuci(2,tt),zonedpctdevse(tt)],0);
    
    [unzonedpctdevmu(tt),unzonedpctdevsigma(tt),unzonedpctdevmuci(:,tt),unzonedpctdevsigmaci(:,tt)]=normfit(unzonedpctdev(tt,:));
    unzonedpctdevse(tt)=unzonedpctdevsigma(tt)/sqrt(length(h));
    unzonedpctdevstats(tt,:)=max([unzonedpctdevmu(tt),unzonedpctdevsigma(tt),unzonedpctdevmuci(1,tt),unzonedpctdevmuci(2,tt),unzonedpctdevse(tt)],0);
    
    [utilmu(tt),utilsigma(tt),utilmuci(:,tt),utilsigmaci(:,tt)]=normfit(pctutilrlts(tt,:));
    utilse(tt)=utilsigma(tt)/sqrt(length(h));
    utilstats(tt,:)=max([utilmu(tt),utilsigma(tt),utilmuci(1,tt),utilmuci(2,tt),utilse(tt)],0);
    
    [devmu(tt),devsigma(tt),devmuci(:,tt),devsigmaci(:,tt)]=normfit(pctdevrlts(tt,:));
    devse(tt)=devsigma(tt)/sqrt(length(h));
    devstats(tt,:)=max([devmu(tt),devsigma(tt),devmuci(1,tt),devmuci(2,tt),devse(tt)],0);
    
    [biglotsmu(tt),biglotssigma(tt),biglotsmuci(:,tt),biglotssigmaci(:,tt)]=normfit(pctbiglotsrlts(tt,:));
    biglotsse(tt)=biglotssigma(tt)/sqrt(length(h));
    biglotsstats(tt,:)=max([biglotsmu(tt),biglotssigma(tt),biglotsmuci(1,tt),biglotsmuci(2,tt),biglotsse(tt)],0);
    
    subagretsold(1:length(h))=agretrlts(1,tt,:);
    [agretsoldmu(tt),agretsoldsigma(tt),agretsoldmuci(:,tt),agretsoldsigmaci(:,tt)]=normfit(subagretsold);
    agretsoldse(tt)=agretsoldsigma(tt)/sqrt(length(h));
    agretsoldstats(tt,:)=max([agretsoldmu(tt),agretsoldsigma(tt),agretsoldmuci(1,tt),agretsoldmuci(2,tt),agretsoldse(tt)],0);
    
    subagretfarm(1:length(h))=agretrlts(2,tt,:);
    [agretfarmmu(tt),agretfarmsigma(tt),agretfarmmuci(:,tt),agretfarmsigmaci(:,tt)]=normfit(subagretfarm);
    agretfarmse(tt)=agretfarmsigma(tt)/sqrt(length(h));
    agretfarmstats(tt,:)=max([agretfarmmu(tt),agretfarmsigma(tt),agretfarmmuci(1,tt),agretfarmmuci(2,tt),agretfarmse(tt)],0);
    
end

zonedenstats=max([zonedenmu,zonedensigma,zonedenmuci(1,:)',zonedenmuci(2,:)',zonedense],0);
[dispmu,dispsigma,dispmuci,dispsigmaci]=normfit(disprlts);
dispse=dispsigma/sqrt(length(h));
dispstats=max([dispmu,dispsigma,dispmuci(1),dispmuci(2),dispse],0);
[distmu,distsigma,distmuci,distsigmaci]=normfit(maxdistrlts);
distse=distsigma/sqrt(length(h));
diststats=max([distmu,distsigma,distmuci(1),distmuci(2),distse],0);

[unzonedpricemu,unzonedpricesigma,unzonedpricemuci,unzonedpricesigmaci]=normfit(zoneavgprice(1,:));
unzonedpricese=unzonedpricesigma/sqrt(length(h));
zoneavgpricestats(1,:)=max([unzonedpricemu,unzonedpricesigma,unzonedpricemuci(1),unzonedpricemuci(2),unzonedpricese],0);
[zonedpricemu,zonedpricesigma,zonedpricemuci,zonedpricesigmaci]=normfit(zoneavgprice(2,:));
zonedpricese=zonedpricesigma/sqrt(length(h));
zoneavgpricestats(2,:)=max([zonedpricemu,zonedpricesigma,zonedpricemuci(1),zonedpricemuci(2),zonedpricese],0);
    
[overplandmu,overplandsigma,overplandmuci,overplandsigmaci]=normfit(overpland(1,:));
overplandse=overplandsigma/sqrt(length(h));
overplandstats(1,:)=max([overplandmu,overplandsigma,overplandmuci(1),overplandmuci(2),overplandse],0);
%     [psoldmu,psoldsigma,psoldmuci,psoldsigmaci]=normfit(plandrlts(1,:));
%     psoldse=psoldsigma/sqrt(length(h));
%     psoldstats=max([psoldmu,psoldsigma,psoldmuci(1),psoldmuci(2),psoldse],0);
%     [pfarmmu,pfarmsigma,pfarmmuci,pfarmsigmaci]=normfit(plandrlts(2,:));
%     pfarmse=pfarmsigma/sqrt(length(h));
%     pfarmstats=max([pfarmmu,pfarmsigma,pfarmmuci(1),pfarmmuci(2),pfarmse],0);



%% Write Excel File Labels %%

rentlabel={'Avg Rents'};
retlabel={'Avg Returns'};
htutillabel={'Avg Utility'};
novaclabel={'No Vacancies'};
htincomelabel={'Income per ht'};
bidlabel={'Share of Bids'};
statlabel={'Mean','Sigma','Low 95 CI','Up 95 CI','SE'};
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
incomelabel={'Mean Income of Occupants'};
northincomelabel={'Mean Income of Occupants in Unzoned Region'};
southincomelabel={'Mean Income of Occupants in Zoned Region'};
outincomelabel={'Mean Income of Ex-occupants'};
zonedpctdevlabel={'Pct. Dev. in Zoned Region'};
unzonedpctdevlabel={'Pct. Dev. in Unzoned Region'};
maxdistlabel={'Max Linear Distance of Development'};
zonepricelabel={'Mean Land Price by Zoning Region'};
northlabel={'North'};
southlabel={'South'};
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
plandtimelabel={'Time Step of Land Sale vs. Land Price ($/ac)'};
tlabel={'Time Step'};
dlabel={'Distance'};
pllabel={'Land Price'};
overplandlabel={'Overall Mean Land Price'};
meanlabel={'Mean'};
sigmalabel={'Sigma'};
offerlabel={'Avg House Offering'};
ideallabel={'Utility Max House Offering'};
proflabel={'Profit Max House Offering'};
difflabel={'Difference'};
varlabel={'Average Broker Variance'};
carrycostlabel={'Carrying Cost'};

%% Write time step file %%
resultsfile=('Results_CHALMS_Coast_0v3.xlsx');
xlswrite(resultsfile,timelabel,'Time Step Stats','B1');
xlswrite(resultsfile,testtime(1),'Time Step Stats','C1');
xlswrite(resultsfile,testtime(2),'Time Step Stats','I1');
xlswrite(resultsfile,testtime(3),'Time Step Stats','O1');
xlswrite(resultsfile,testtime(4),'Time Step Stats','U1');
xlswrite(resultsfile,testtime(5),'Time Step Stats','AA1');
xlswrite(resultsfile,lottypelabel,'Time Step Stats','B2');
xlswrite(resultsfile,statlabel,'Time Step Stats','C2');
xlswrite(resultsfile,statlabel,'Time Step Stats','I2');
xlswrite(resultsfile,statlabel,'Time Step Stats','O2');
xlswrite(resultsfile,statlabel,'Time Step Stats','U2');
xlswrite(resultsfile,statlabel,'Time Step Stats','AA2');
xlswrite(resultsfile,rentlabel,'Time Step Stats','A3');
xlswrite(resultsfile,lotlabel,'Time Step Stats','A22');
xlswrite(resultsfile,lottypeset,'Time Step Stats','B3');
xlswrite(resultsfile,lottypeset,'Time Step Stats','B22');
xlswrite(resultsfile,lotsumlabel,'Time Step Stats','B40');
xlswrite(resultsfile,statlabel,'Time Step Stats','C42');
xlswrite(resultsfile,distlabel,'Time Step Stats','B42');
xlswrite(resultsfile,zonedenlabel,'Time Step Stats','A42');
xlswrite(resultsfile,plandtimelabel,'Time Step Stats','B51');
xlswrite(resultsfile,tlabel,'Time Step Stats','A52');
xlswrite(resultsfile,pllabel,'Time Step Stats','A53');
xlswrite(resultsfile,planddistlabel,'Time Step Stats','B55');
xlswrite(resultsfile,dlabel,'Time Step Stats','A56');
xlswrite(resultsfile,pllabel,'Time Step Stats','A57');

xlswrite(resultsfile,rentstats(:,:,1),'Time Step Stats','C3');
xlswrite(resultsfile,rentstats(:,:,2),'Time Step Stats','I3');
xlswrite(resultsfile,rentstats(:,:,3),'Time Step Stats','O3');
xlswrite(resultsfile,rentstats(:,:,4),'Time Step Stats','U3');
xlswrite(resultsfile,rentstats(:,:,5),'Time Step Stats','AA3');
xlswrite(resultsfile,lotstats(:,:,1),'Time Step Stats','C22');
xlswrite(resultsfile,lotstats(:,:,2),'Time Step Stats','I22');
xlswrite(resultsfile,lotstats(:,:,3),'Time Step Stats','O22');
xlswrite(resultsfile,lotstats(:,:,4),'Time Step Stats','U22');
xlswrite(resultsfile,lotstats(:,:,5),'Time Step Stats','AA22');
xlswrite(resultsfile,lotsumstats(1,:),'Time Step Stats','C40');
xlswrite(resultsfile,lotsumstats(2,:),'Time Step Stats','I40');
xlswrite(resultsfile,lotsumstats(3,:),'Time Step Stats','O40');
xlswrite(resultsfile,lotsumstats(4,:),'Time Step Stats','U40');
xlswrite(resultsfile,lotsumstats(5,:),'Time Step Stats','AA40');

xlswrite(resultsfile,zonedenstats,'Time Step Stats','C43');
xlswrite(resultsfile,zonedist,'Time Step Stats','B43');
xlswrite(resultsfile,maxdistlabel,'Time Step Stats','I42');
xlswrite(resultsfile,statlabel,'Time Step Stats','I43');
xlswrite(resultsfile,diststats,'Time Step Stats','I44');
xlswrite(resultsfile,displabel,'Time Step Stats','I46');
xlswrite(resultsfile,statlabel,'Time Step Stats','I47');
xlswrite(resultsfile,dispstats,'Time Step Stats','I48');

xlswrite(resultsfile,avgplandtime','Time Step Stats','B52');
xlswrite(resultsfile,avgplanddist','Time Step Stats','B56');
xlswrite(resultsfile,carrycostlabel,'Time Step Stats','A54');
xlswrite(resultsfile,carrycostrlts(1,TSTART+1:TMAX),'Time Step Stats','B54');

xlswrite(resultsfile,testtime','Landscape Results','A3');
xlswrite(resultsfile,devlabel,'Landscape Results','B1');
xlswrite(resultsfile,statlabel,'Landscape Results','B2');
xlswrite(resultsfile,devstats,'Landscape Results','B3');
xlswrite(resultsfile,testtime','Landscape Results','A11');
xlswrite(resultsfile,vaclabel,'Landscape Results','B9');
xlswrite(resultsfile,statlabel,'Landscape Results','B10');
xlswrite(resultsfile,vacstats,'Landscape Results','B11');
% xlswrite(resultsfile,testtime','Landscape Results','A20');
xlswrite(resultsfile,farmunsoldlabel,'Landscape Results','A20');
xlswrite(resultsfile,farmsoldlabel,'Landscape Results','A21');
xlswrite(resultsfile,pfarmlabel,'Landscape Results','B18');
xlswrite(resultsfile,statlabel,'Landscape Results','B19');
xlswrite(resultsfile,agretfarmstats(5,:),'Landscape Results','B20');
% xlswrite(resultsfile,testtime','Landscape Results','A28');
% xlswrite(resultsfile,psoldlabel,'Landscape Results','B26');
% xlswrite(resultsfile,statlabel,'Landscape Results','B27');
xlswrite(resultsfile,agretsoldstats(5,:),'Landscape Results','B21');
xlswrite(resultsfile,testtime','Landscape Results','H3');
xlswrite(resultsfile,unzonedpctdevlabel,'Landscape Results','I1');
xlswrite(resultsfile,statlabel,'Landscape Results','I2');
xlswrite(resultsfile,unzonedpctdevstats,'Landscape Results','I3');
xlswrite(resultsfile,testtime','Landscape Results','H11');
xlswrite(resultsfile,zonedpctdevlabel,'Landscape Results','I9');
xlswrite(resultsfile,statlabel,'Landscape Results','I10');
xlswrite(resultsfile,zonedpctdevstats,'Landscape Results','I11');
xlswrite(resultsfile,testtime','Landscape Results','H20');
xlswrite(resultsfile,incomelabel,'Landscape Results','I18');
xlswrite(resultsfile,statlabel,'Landscape Results','I19');
xlswrite(resultsfile,incomestats,'Landscape Results','I20');
xlswrite(resultsfile,testtime','Landscape Results','H28');
xlswrite(resultsfile,outincomelabel,'Landscape Results','I26');
xlswrite(resultsfile,statlabel,'Landscape Results','I27');
xlswrite(resultsfile,outincomestats,'Landscape Results','I28');
xlswrite(resultsfile,zonepricelabel,'Landscape Results','B24');
xlswrite(resultsfile,statlabel,'Landscape Results','B25');
xlswrite(resultsfile,northlabel,'Landscape Results','A26');
xlswrite(resultsfile,southlabel,'Landscape Results','A27');
xlswrite(resultsfile,zoneavgpricestats,'Landscape Results','B26');
xlswrite(resultsfile,overplandlabel,'Landscape Results','B29');
xlswrite(resultsfile,statlabel,'Landscape Results','B30');
xlswrite(resultsfile,overplandstats,'Landscape Results','B31');
% xlswrite(resultsfile,zonepricelabel,'Landscape Results','B17');
% xlswrite(resultsfile,statlabel,'Landscape Results','B18');
% xlswrite(resultsfile,northlabel,'Landscape Results','A19');
% xlswrite(resultsfile,southlabel,'Landscape Results','A20');
% xlswrite(resultsfile,zoneavgpricestats,'Landscape Results','B19');
xlswrite(resultsfile,northincomelabel,'Landscape Results','P1');
xlswrite(resultsfile,statlabel,'Landscape Results','P2');
xlswrite(resultsfile,testtime','Landscape Results','O3');
xlswrite(resultsfile,northincomestats,'Landscape Results','P3');
xlswrite(resultsfile,southincomelabel,'Landscape Results','P9');
xlswrite(resultsfile,statlabel,'Landscape Results','P10');
xlswrite(resultsfile,testtime','Landscape Results','O11');
xlswrite(resultsfile,southincomestats,'Landscape Results','P11');
xlswrite(resultsfile,distlabel,'Landscape Results','B34');
xlswrite(resultsfile,lotsizelabel,'Landscape Results','A36');
xlswrite(resultsfile,finedistances,'Landscape Results','B35');
xlswrite(resultsfile,meanlotsize,'Landscape Results','B36');
xlswrite(resultsfile,hpaclabel,'Landscape Results','A37');
xlswrite(resultsfile,meanhouseperacre,'Landscape Results','B37');
xlswrite(resultsfile,pllabel,'Landscape Results','A38');
xlswrite(resultsfile,meanplandint,'Landscape Results','B38');


xlswrite(resultsfile,timelabel,'North Results','B1');
xlswrite(resultsfile,testtime(1),'North Results','C1');
xlswrite(resultsfile,testtime(2),'North Results','I1');
xlswrite(resultsfile,testtime(3),'North Results','O1');
xlswrite(resultsfile,testtime(4),'North Results','U1');
xlswrite(resultsfile,testtime(5),'North Results','AA1');
xlswrite(resultsfile,lottypelabel,'North Results','B2');
xlswrite(resultsfile,statlabel,'North Results','C2');
xlswrite(resultsfile,statlabel,'North Results','I2');
xlswrite(resultsfile,statlabel,'North Results','O2');
xlswrite(resultsfile,statlabel,'North Results','U2');
xlswrite(resultsfile,statlabel,'North Results','AA2');
xlswrite(resultsfile,rentlabel,'North Results','A3');
xlswrite(resultsfile,lotlabel,'North Results','A22');
xlswrite(resultsfile,lottypeset,'North Results','B3');
xlswrite(resultsfile,lottypeset,'North Results','B22');
xlswrite(resultsfile,nrentstats(:,:,1),'North Results','C3');
xlswrite(resultsfile,nrentstats(:,:,2),'North Results','I3');
xlswrite(resultsfile,nrentstats(:,:,3),'North Results','O3');
xlswrite(resultsfile,nrentstats(:,:,4),'North Results','U3');
xlswrite(resultsfile,nrentstats(:,:,5),'North Results','AA3');
xlswrite(resultsfile,unzonedlotstats(:,:,1),'North Results','C22');
xlswrite(resultsfile,unzonedlotstats(:,:,2),'North Results','I22');
xlswrite(resultsfile,unzonedlotstats(:,:,3),'North Results','O22');
xlswrite(resultsfile,unzonedlotstats(:,:,4),'North Results','U22');
xlswrite(resultsfile,unzonedlotstats(:,:,5),'North Results','AA22');

xlswrite(resultsfile,timelabel,'South Results','B1');
xlswrite(resultsfile,testtime(1),'South Results','C1');
xlswrite(resultsfile,testtime(2),'South Results','I1');
xlswrite(resultsfile,testtime(3),'South Results','O1');
xlswrite(resultsfile,testtime(4),'South Results','U1');
xlswrite(resultsfile,testtime(5),'South Results','AA1');
xlswrite(resultsfile,lottypelabel,'South Results','B2');
xlswrite(resultsfile,statlabel,'South Results','C2');
xlswrite(resultsfile,statlabel,'South Results','I2');
xlswrite(resultsfile,statlabel,'South Results','O2');
xlswrite(resultsfile,statlabel,'South Results','U2');
xlswrite(resultsfile,statlabel,'South Results','AA2');
xlswrite(resultsfile,rentlabel,'South Results','A3');
xlswrite(resultsfile,lotlabel,'South Results','A22');
xlswrite(resultsfile,lottypeset,'South Results','B3');
xlswrite(resultsfile,lottypeset,'South Results','B22');
xlswrite(resultsfile,srentstats(:,:,1),'South Results','C3');
xlswrite(resultsfile,srentstats(:,:,2),'South Results','I3');
xlswrite(resultsfile,srentstats(:,:,3),'South Results','O3');
xlswrite(resultsfile,srentstats(:,:,4),'South Results','U3');
xlswrite(resultsfile,srentstats(:,:,5),'South Results','AA3');
xlswrite(resultsfile,zonedlotstats(:,:,1),'South Results','C22');
xlswrite(resultsfile,zonedlotstats(:,:,2),'South Results','I22');
xlswrite(resultsfile,zonedlotstats(:,:,3),'South Results','O22');
xlswrite(resultsfile,zonedlotstats(:,:,4),'South Results','U22');
xlswrite(resultsfile,zonedlotstats(:,:,5),'South Results','AA22');

obslabel={'Obs.'};
avgmaplabel={'Avg Map'};
pctdevlabel={'%dev'};
threshlabel={'thresh'};
xlswrite(resultsfile,obslabel,'AvgMap15','B25');
xlswrite(resultsfile,avgmaplabel,'AvgMap15','C25');
xlswrite(resultsfile,pctdevlabel,'AvgMap15','A26');
xlswrite(resultsfile,threshlabel,'AvgMap15','A27');
xlswrite(resultsfile,devstats(2,1),'AvgMap15','B26');
xlswrite(resultsfile,testpctdev(2),'AvgMap15','C26');
xlswrite(resultsfile,avgthresh(2),'AvgMap15','B27');

xlswrite(resultsfile,obslabel,'AvgMap20','B25');
xlswrite(resultsfile,avgmaplabel,'AvgMap20','C25');
xlswrite(resultsfile,pctdevlabel,'AvgMap20','A26');
xlswrite(resultsfile,threshlabel,'AvgMap20','A27');
xlswrite(resultsfile,devstats(3,1),'AvgMap20','B26');
xlswrite(resultsfile,testpctdev(3),'AvgMap20','C26');
xlswrite(resultsfile,avgthresh(3),'AvgMap20','B27');

xlswrite(resultsfile,obslabel,'AvgMap25','B25');
xlswrite(resultsfile,avgmaplabel,'AvgMap25','C25');
xlswrite(resultsfile,pctdevlabel,'AvgMap25','A26');
xlswrite(resultsfile,threshlabel,'AvgMap25','A27');
xlswrite(resultsfile,devstats(4,1),'AvgMap25','B26');
xlswrite(resultsfile,testpctdev(4),'AvgMap25','C26');
xlswrite(resultsfile,avgthresh(4),'AvgMap25','B27');

xlswrite(resultsfile,obslabel,'AvgMap30','B25');
xlswrite(resultsfile,avgmaplabel,'AvgMap30','C25');
xlswrite(resultsfile,pctdevlabel,'AvgMap30','A26');
xlswrite(resultsfile,threshlabel,'AvgMap30','A27');
xlswrite(resultsfile,devstats(5,1),'AvgMap30','B26');
xlswrite(resultsfile,testpctdev(5),'AvgMap30','C26');
xlswrite(resultsfile,avgthresh(5),'AvgMap30','B27');

xlswrite(resultsfile,timelabel,'Return Stats','B1');
xlswrite(resultsfile,(TSTART+1:TMAX),'Return Stats','C1');
xlswrite(resultsfile,lottypelabel,'Return Stats','B2');
xlswrite(resultsfile,retlabel,'Return Stats','A3');
xlswrite(resultsfile,sigmalabel,'Return Stats','A22');
xlswrite(resultsfile,lottypeset,'Return Stats','B3');
xlswrite(resultsfile,lottypeset,'Return Stats','B22');
xlswrite(resultsfile,meanlabel,'Return Stats','C2');
xlswrite(resultsfile,reshape(retstats(:,1,:),[HT 20]),'Return Stats','C3');
xlswrite(resultsfile,reshape(retstats(:,2,:),[HT 20]),'Return Stats','C22');

% xlswrite(resultsfile,timelabel,'HtUtility Stats','B1');
% xlswrite(resultsfile,(TSTART+1:TMAX),'HtUtility Stats','C1');
% xlswrite(resultsfile,lottypelabel,'HtUtility Stats','B2');
% xlswrite(resultsfile,htutillabel,'HtUtility Stats','A3');
% xlswrite(resultsfile,novaclabel,'HtUtility Stats','A22');
% xlswrite(resultsfile,lottypeset,'HtUtility Stats','B3');
% xlswrite(resultsfile,lottypeset,'HtUtility Stats','B22');
% xlswrite(resultsfile,meanlabel,'HtUtility Stats','C2');
% xlswrite(resultsfile,reshape(htutilvacstats(:,1,:),[HT 20]),'HtUtility Stats','C3');
% xlswrite(resultsfile,reshape(htutilnovacstats(:,1,:),[HT 20]),'HtUtility Stats','C22');

xlswrite(resultsfile,timelabel,'HtIncome Stats','B1');
xlswrite(resultsfile,(TSTART+1:TMAX),'HtIncome Stats','C1');
xlswrite(resultsfile,lottypelabel,'HtIncome Stats','B2');
xlswrite(resultsfile,htincomelabel,'HtIncome Stats','A3');
xlswrite(resultsfile,sigmalabel,'HtIncome Stats','A22');
xlswrite(resultsfile,lottypeset,'HtIncome Stats','B3');
xlswrite(resultsfile,lottypeset,'HtIncome Stats','B22');
xlswrite(resultsfile,meanlabel,'HtIncome Stats','C2');
xlswrite(resultsfile,reshape(htincomestats(:,1,:),[HT 20]),'HtIncome Stats','C3');
xlswrite(resultsfile,reshape(htincomestats(:,2,:),[HT 20]),'HtIncome Stats','C22');

xlswrite(resultsfile,timelabel,'BidShare Stats','B1');
xlswrite(resultsfile,(TSTART+1:TMAX),'BidShare Stats','C1');
xlswrite(resultsfile,lottypelabel,'BidShare Stats','B2');
xlswrite(resultsfile,bidlabel,'BidShare Stats','A3');
xlswrite(resultsfile,sigmalabel,'BidShare Stats','A22');
xlswrite(resultsfile,lottypeset,'BidShare Stats','B3');
xlswrite(resultsfile,lottypeset,'BidShare Stats','B22');
xlswrite(resultsfile,meanlabel,'BidShare Stats','C2');
xlswrite(resultsfile,reshape(bidsharestats(:,1,:),[HT 20]),'BidShare Stats','C3');
xlswrite(resultsfile,reshape(bidsharestats(:,2,:),[HT 20]),'BidShare Stats','C22');

xlswrite(resultsfile,offerlabel,'HouseOfferings','A2');
xlswrite(resultsfile,timelabel,'HouseOfferings','D3');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','A4');
xlswrite(resultsfile,lottypeset,'HouseOfferings','A5');
xlswrite(resultsfile,testtime,'HouseOfferings','B4');
xlswrite(resultsfile,offerset,'HouseOfferings','B5');

xlswrite(resultsfile,varlabel,'HouseOfferings','H2');
xlswrite(resultsfile,timelabel,'HouseOfferings','K3');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','H4');
xlswrite(resultsfile,lottypeset,'HouseOfferings','H5');
xlswrite(resultsfile,testtime,'HouseOfferings','I4');
xlswrite(resultsfile,avgvar(:,testtime-9),'HouseOfferings','I5');

xlswrite(resultsfile,ideallabel,'HouseOfferings','A24');
xlswrite(resultsfile,timelabel,'HouseOfferings','D25');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','A26');
xlswrite(resultsfile,lottypeset,'HouseOfferings','A27');
xlswrite(resultsfile,testtime,'HouseOfferings','B26');
xlswrite(resultsfile,avgidealset,'HouseOfferings','B27');

xlswrite(resultsfile,proflabel,'HouseOfferings','A46');
xlswrite(resultsfile,timelabel,'HouseOfferings','D47');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','A48');
xlswrite(resultsfile,lottypeset,'HouseOfferings','A49');
xlswrite(resultsfile,testtime,'HouseOfferings','B48');
xlswrite(resultsfile,avgprofset,'HouseOfferings','B49');

xlswrite(resultsfile,difflabel,'HouseOfferings','H24');
xlswrite(resultsfile,timelabel,'HouseOfferings','K25');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','H26');
xlswrite(resultsfile,lottypeset,'HouseOfferings','H27');
xlswrite(resultsfile,testtime,'HouseOfferings','I26');
xlswrite(resultsfile,diffset,'HouseOfferings','I27');

xlswrite(resultsfile,difflabel,'HouseOfferings','H46');
xlswrite(resultsfile,timelabel,'HouseOfferings','K47');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','H48');
xlswrite(resultsfile,lottypeset,'HouseOfferings','H49');
xlswrite(resultsfile,testtime,'HouseOfferings','I48');
xlswrite(resultsfile,diffprofset,'HouseOfferings','I49');

