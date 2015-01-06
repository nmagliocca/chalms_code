%%%%%%%% Write results to Excel file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MRUNS=30;
TMAX=30;
TSTART=10;
NPARMS=1;
EXPTRUNS=4;
ERUNS=EXPTRUNS^NPARMS;
batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
    repmat((1:MRUNS)',ERUNS,1)];
batchruns=mat2cell(reshape(1:MRUNS*ERUNS,MRUNS,ERUNS),MRUNS,ones(1,ERUNS));

NWIDTH=80;
NLENGTH=80;
NCELLS=NLENGTH*NWIDTH;
HT=8;
discount=0.05;

VARLAYER=zeros(NLENGTH*NWIDTH,EXPTRUNS*MRUNS);

BASELAYER=zeros(NCELLS,1);
COAST=zeros(NLENGTH,NWIDTH);
COAST(:,1:5)=1;
icoast=find(COAST==1);
SCAPE=ones(NLENGTH,NWIDTH);
iscape=(SCAPE==1);
iscapelist=find(iscape==1);
% COAST=zeros(NLENGTH,NWIDTH);
% SCAPE=zeros(NLENGTH,NWIDTH);
% COAST(:,1:5)=1;
% SCAPE(COAST~=1)=1;
% DISTANCE=ones(NLENGTH,NWIDTH);
% dist2hznnei=zeros(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
% dist2vrtnei=zeros(NLENGTH,NWIDTH);
NZONES=25;
ZONES=cell(NZONES,3);   %[cellid,min zone,max zone]
% ZONES=zeros(NLENGTH,NWIDTH);
zoning=zeros(NZONES,2);     %[(min lotsize) (max lotsize)]
ZONEMAP=zeros(NLENGTH,NWIDTH);
CZONEMAP=zeros(NLENGTH,NWIDTH); %coatal zone = elevated expected land values
nzoneslong=sqrt(NZONES);
nzoneswide=sqrt(NZONES);
extralengthz=rem(NLENGTH,nzoneslong);
extrawidthz=rem(NWIDTH,nzoneswide);
izonelength=(NLENGTH-extralengthz)/nzoneslong;
izonewidth=(NWIDTH-extrawidthz)/nzoneswide;
izoneborder=zeros(1,NWIDTH);
borderdist=zeros(1,NWIDTH);
% nzoneslong=sqrt(NZONES);
% nzoneswide=sqrt(NZONES);
% extralengthz=rem(NLENGTH,nzoneslong);
% extrawidthz=rem(NWIDTH,nzoneswide);
% izonelength=(length(BASELAYER(:,1))-extralengthz)/nzoneslong;
% izonewidth=(length(BASELAYER(1,:))-extrawidthz)/nzoneswide;

zonemarklong=1;
zonemarkwide=1;
for ii=1:nzoneswide
    for jj=1:nzoneslong
        Zones(zonemarklong:zonemarklong+izonelength-1,...
            zonemarkwide:zonemarkwide+izonewidth-1)=...
            ii*nzoneslong-(nzoneslong-jj);        
        
        if jj==nzoneslong && extralengthz > 0
            Zones(izonelength*jj+1:izonelength*jj+extralengthz,...
                zonemarkwide:zonemarkwide+izonewidth-1)=...
                ii*nzoneslong-(nzoneslong-jj);        
        end
        if ii==nzoneswide && extrawidthz > 0
            Zones(zonemarklong:zonemarklong+izonelength-1,...
                izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
                ii*nzoneslong-(nzoneslong-jj);
        end
        zonemarklong=mod(izonelength*jj+1,izonelength*nzoneslong);
    end     
    if jj==nzoneslong && ii==nzoneswide && extralengthz > 0
        Zones(izonelength*jj+1:izonelength*jj+extralengthz,...
            izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
            ii*nzoneslong-(nzoneslong-jj);
    end
    zonemarkwide=mod(izonewidth*ii+1,izonewidth*nzoneswide);
end
for iz=1:NZONES
    ZONES{iz,1}=find(Zones==iz);
end
izoned=cat(1,ZONES{[3:5 9:10 14:15 19:20 23:25],1});
inotzoned=cat(1,ZONES{[1:2 6:8 11:13 16:18 21:22],1});
% izoned=ismember(ZONES,[3:5 9:10 14:15 19:20 23:25]);
% inotzoned=~ismember(ZONES,[3:5 9:10 14:15 19:20 23:25]);
ZONEMAP(izoned)=1;
CZONEMAP(:,1:50)=1;
CZONEMAP(:,51:70)=2;
CZONEMAP(:,71:NWIDTH)=3;
unzonedacres=length(find(ZONEMAP==0));
zonedacres=length(find(ZONEMAP==1));
inlandacres=length(find(CZONEMAP==1));
middleacres=length(find(CZONEMAP==2));
coastalacres=length(find(CZONEMAP==3));
icoastal=find(CZONEMAP==3);
imiddle=find(CZONEMAP==2);
iinland=find(CZONEMAP==1);

northindex=find(ZONEMAP==0);
southindex=find(ZONEMAP==1);

cd X:\model_results\CHALMS_alt_storm_climate
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('storm_clim3',fnamescell(1,:),11);
hind=find(h==1);
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files
load DIST2CBD_east
cd X:\model_results\CHALMS_alt_storm_climate
distpt=ceil(min(min(dist2cbd)):max(max(dist2cbd)));

for i=1:NWIDTH
    izoneborder(i)=find(ZONEMAP(:,i)==1,1,'first');
    borderdist(i)=dist2cbd(izoneborder(i),i)*0.0395;
end
coastdist=cumsum(reshape(SCAPE,NLENGTH,NWIDTH),2);
coastprox=max(max(coastdist))-coastdist;
%%%%%
MILEDIST=dist2cbd.*0.0395;
CDIST=coastprox*0.0395;
%%%%%

testtime=[10 15 20 25 30];
zonedist=[0 14 28 42 56 70 84]'.*0.0395;
timeset=(TSTART+1:TMAX);
z = 1000*[0.00025    1.5000
    0.00025    2.5000
    0.0005    1.5000
    0.0005    2.5000
    0.0010    1.5000
    0.0010    2.5000
    0.0020    1.5000
    0.0020    2.5000];

lotsizemap=zeros(NLENGTH,NWIDTH);

% fnames=dir;
% fnamescell=struct2cell(fnames);
%  %%% CHANGE FILE NAME %%%%
% hind=strmatch('results_CHALMS_Coast_base',fnamescell(1,:));



bidsharerlts=zeros(HT,length(TSTART+1:TMAX),length(hind));
bidsharestats=zeros(HT,5,length(TSTART+1:TMAX));
bidsharemu=zeros(HT,length(TSTART+1:TMAX));
bidsharesigma=zeros(HT,length(TSTART+1:TMAX));
bidsharemuci=zeros(2,HT,length(TSTART+1:TMAX));
bidsharesigmaci=zeros(2,HT,length(TSTART+1:TMAX));
bidsharese=zeros(HT,length(TSTART+1:TMAX));

htutilvacrlts=zeros(HT,length(TSTART+1:TMAX),length(hind));
htutilvacstats=zeros(HT,5,length(TSTART+1:TMAX));
htutilvacmu=zeros(HT,length(TSTART+1:TMAX));
htutilvacsigma=zeros(HT,length(TSTART+1:TMAX));
htutilvacmuci=zeros(2,HT,length(TSTART+1:TMAX));
htutilvacsigmaci=zeros(2,HT,length(TSTART+1:TMAX));
htutilvacse=zeros(HT,length(TSTART+1:TMAX));

htutilnovacrlts=zeros(HT,length(TSTART+1:TMAX),length(hind));
htutilnovacstats=zeros(HT,5,length(TSTART+1:TMAX));
htutilnovacmu=zeros(HT,length(TSTART+1:TMAX));
htutilnovacsigma=zeros(HT,length(TSTART+1:TMAX));
htutilnovacmuci=zeros(2,HT,length(TSTART+1:TMAX));
htutilnovacsigmaci=zeros(2,HT,length(TSTART+1:TMAX));
htutilnovacse=zeros(HT,length(TSTART+1:TMAX));

htincomerlts=zeros(HT,length(TSTART+1:TMAX),length(hind));
htincomestats=zeros(HT,5,length(TSTART+1:TMAX));
htincomemu=zeros(HT,length(TSTART+1:TMAX));
htincomesigma=zeros(HT,length(TSTART+1:TMAX));
htincomemuci=zeros(2,HT,length(TSTART+1:TMAX));
htincomesigmaci=zeros(2,HT,length(TSTART+1:TMAX));
htincomese=zeros(HT,length(TSTART+1:TMAX));

retrlts=zeros(HT,length(TSTART+1:TMAX),length(hind));
retstats=zeros(HT,5,length(TSTART+1:TMAX));
retmu=zeros(HT,length(TSTART+1:TMAX));
retsigma=zeros(HT,length(TSTART+1:TMAX));
retmuci=zeros(2,HT,length(TSTART+1:TMAX));
retsigmaci=zeros(2,HT,length(TSTART+1:TMAX));
retse=zeros(HT,length(TSTART+1:TMAX));

fullrentrlts=zeros(HT,TMAX,length(hind));
rentrlts=zeros(HT,length(testtime),length(hind));
rentstats=zeros(HT,5,length(testtime));
rentmu=zeros(HT,length(testtime));
rentsigma=zeros(HT,length(testtime));
rentmuci=zeros(2,HT,length(testtime));
rentsigmaci=zeros(2,HT,length(testtime));
rentse=zeros(HT,length(testtime));


zonepctdev=zeros(2,length(testtime));
zonenumlots=zeros(2,length(testtime));
% zoneavgprice=zeros(2,length(testtime));
czonepctdev=zeros(2,length(testtime));
czonenumlots=zeros(2,length(testtime));
% czoneavgprice=zeros(2,length(testtime));
% incomestats=zeros(1,length(testtime));

northrents=zeros(HT,length(testtime),length(hind));
nrentstats=zeros(HT,5,length(testtime));
nrentmu=zeros(HT,length(testtime));
nrentsigma=zeros(HT,length(testtime));
nrentmuci=zeros(2,HT,length(testtime));
nrentsigmaci=zeros(2,HT,length(testtime));
nrentse=zeros(HT,length(testtime));

southrents=zeros(HT,length(testtime),length(hind));
srentstats=zeros(HT,5,length(testtime));
srentmu=zeros(HT,length(testtime));
srentsigma=zeros(HT,length(testtime));
srentmuci=zeros(2,HT,length(testtime));
srentsigmaci=zeros(2,HT,length(testtime));
srentse=zeros(HT,length(testtime));

coastrents=zeros(HT,length(testtime),length(hind));
crentstats=zeros(HT,5,length(testtime));
crentmu=zeros(HT,length(testtime));
crentsigma=zeros(HT,length(testtime));
crentmuci=zeros(2,HT,length(testtime));
crentsigmaci=zeros(2,HT,length(testtime));
crentse=zeros(HT,length(testtime));

middlerents=zeros(HT,length(testtime),length(hind));
midrentstats=zeros(HT,5,length(testtime));
midrentmu=zeros(HT,length(testtime));
midrentsigma=zeros(HT,length(testtime));
midrentmuci=zeros(2,HT,length(testtime));
midrentsigmaci=zeros(2,HT,length(testtime));
midrentse=zeros(HT,length(testtime));

inlandrents=zeros(HT,length(testtime),length(hind));
inrentstats=zeros(HT,5,length(testtime));
inrentmu=zeros(HT,length(testtime));
inrentsigma=zeros(HT,length(testtime));
inrentmuci=zeros(2,HT,length(testtime));
inrentsigmaci=zeros(2,HT,length(testtime));
inrentse=zeros(HT,length(testtime));

vacrlts=zeros(length(testtime),length(hind));
pctutilrlts=zeros(length(testtime),length(hind));
numlotrlts=zeros(HT,length(testtime),length(hind));
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

coastlotstats=zeros(HT,5,length(testtime));
coastlotmu=zeros(HT,length(testtime));
coastlotsigma=zeros(HT,length(testtime));
coastlotmuci=zeros(2,HT,length(testtime));
coastlotsigmaci=zeros(2,HT,length(testtime));
coastlotse=zeros(HT,length(testtime));

middlelotstats=zeros(HT,5,length(testtime));
middlelotmu=zeros(HT,length(testtime));
middlelotsigma=zeros(HT,length(testtime));
middlelotmuci=zeros(2,HT,length(testtime));
middlelotsigmaci=zeros(2,HT,length(testtime));
middlelotse=zeros(HT,length(testtime));

inlandlotstats=zeros(HT,5,length(testtime));
inlandlotmu=zeros(HT,length(testtime));
inlandlotsigma=zeros(HT,length(testtime));
inlandlotmuci=zeros(2,HT,length(testtime));
inlandlotsigmaci=zeros(2,HT,length(testtime));
inlandlotse=zeros(HT,length(testtime));

lotsumrlts=zeros(length(testtime),length(hind));
disprlts=zeros(1,length(hind));
maxdistrlts=zeros(1,length(hind));
pctdevrlts=zeros(length(testtime),length(hind));
testpctdev=zeros(length(testtime),1);
pctbiglotsrlts=zeros(length(testtime),length(hind));
plandrlts=zeros(50,4,length(hind));
plandinfo=zeros([],5); %[farmerid selltime sellprice dist acres]
avgplandtime=zeros(20,2);
avgplanddist=zeros(50,2);
agretrlts=zeros(2,length(testtime),length(hind));
zonedenrlts=zeros(HT,length(hind));
zonedenmu=zeros(HT,1);
zonedensigma=zeros(HT,1);
zonedenmuci=zeros(2,HT);
zonedensigmaci=zeros(2,HT);
zonedense=zeros(HT,1);
czonedenrlts=zeros(HT,length(hind));
czonedenmu=zeros(HT,1);
czonedensigma=zeros(HT,1);
czonedenmuci=zeros(2,HT);
czonedensigmaci=zeros(2,HT);
czonedense=zeros(HT,1);

zoneavgpricestats=zeros(2,5);
czoneavgpricestats=zeros(3,5);
overplandstats=zeros(1,5);

unzonedpctdev=zeros(length(testtime),length(hind));
zonedpctdev=zeros(length(testtime),length(hind));
zonednumlots=zeros(HT,length(testtime),length(hind));
unzonednumlots=zeros(HT,length(testtime),length(hind));
zoneavgprice=zeros(2,length(hind));
coastalpctdev=zeros(length(testtime),length(hind));
middlepctdev=zeros(length(testtime),length(hind));
inlandpctdev=zeros(length(testtime),length(hind));
inlandnumlots=zeros(HT,length(testtime),length(hind));
coastalnumlots=zeros(HT,length(testtime),length(hind));
middlenumlots=zeros(HT,length(testtime),length(hind));
czoneavgprice=zeros(3,length(hind));
overpland=zeros(1,length(hind));
incomerlts=zeros(length(testtime),length(hind));
outincomerlts=zeros(length(testtime),length(hind));

northincome=zeros(length(testtime),length(hind));
southincome=zeros(length(testtime),length(hind));
coastalincome=zeros(length(testtime),length(hind));
middleincome=zeros(length(testtime),length(hind));
inlandincome=zeros(length(testtime),length(hind));

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

coastpctdevstats=zeros(length(testtime),5);
coastpctdevmu=zeros(length(testtime),1);
coastpctdevsigma=zeros(length(testtime),1);
coastpctdevmuci=zeros(2,length(testtime));
coastpctdevsigmaci=zeros(2,length(testtime));
coastpctdevse=zeros(length(testtime),1);

middlepctdevstats=zeros(length(testtime),5);
middlepctdevmu=zeros(length(testtime),1);
middlepctdevsigma=zeros(length(testtime),1);
middlepctdevmuci=zeros(2,length(testtime));
middlepctdevsigmaci=zeros(2,length(testtime));
middlepctdevse=zeros(length(testtime),1);

inlandpctdevstats=zeros(length(testtime),5);
inlandpctdevmu=zeros(length(testtime),1);
inlandpctdevsigma=zeros(length(testtime),1);
inlandpctdevmuci=zeros(2,length(testtime));
inlandpctdevsigmaci=zeros(2,length(testtime));
inlandpctdevse=zeros(length(testtime),1);

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

coastincomestats=zeros(length(testtime),5);
coastincomemu=zeros(length(testtime),1);
coastincomesigma=zeros(length(testtime),1);
coastincomemuci=zeros(2,length(testtime));
coastincomesigmaci=zeros(2,length(testtime));
coastincomese=zeros(length(testtime),1);

middleincomestats=zeros(length(testtime),5);
middleincomemu=zeros(length(testtime),1);
middleincomesigma=zeros(length(testtime),1);
middleincomemuci=zeros(2,length(testtime));
middleincomesigmaci=zeros(2,length(testtime));
middleincomese=zeros(length(testtime),1);

inlandincomestats=zeros(length(testtime),5);
inlandincomemu=zeros(length(testtime),1);
inlandincomesigma=zeros(length(testtime),1);
inlandincomemuci=zeros(2,length(testtime));
inlandincomesigmaci=zeros(2,length(testtime));
inlandincomese=zeros(length(testtime),1);

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

Exptretrlts=zeros(HT,TMAX,length(hind));
Realretrlts=zeros(HT,TMAX,length(hind));
Nlots=zeros(HT,TMAX,length(hind));
IDEALSET=zeros(HT,TMAX,length(hind));
PROFSET=zeros(HT,TMAX,length(hind));
AVGBVAR=zeros(HT,TMAX,length(hind));

lottypemap=zeros(NLENGTH,NWIDTH,HT,length(testtime));
secprobmap=zeros(NLENGTH,NWIDTH,HT,length(testtime));
subvarmap=zeros(NLENGTH,NWIDTH,length(testtime));

VARMAP=zeros(NLENGTH,NWIDTH,length(testtime));

VARLAYERSAVE=zeros(NLENGTH,NWIDTH);
avgthresh=[1 0.36 0.4 0.4 0.45];
% avgthresh=[1 0.33 0.4 0.47 0.63];
AVGMAP=zeros(NLENGTH,NWIDTH,length(testtime));
CONFMAP=zeros(NLENGTH,NWIDTH,length(testtime));
carcost=zeros(TMAX,length(hind));
carrycostrlts=zeros(1,TMAX);

PMAPS=zeros(NLENGTH*NWIDTH,length(testtime),length(hind));
avgPREFMAP=zeros(NLENGTH*NWIDTH,length(testtime));

t30maps=zeros(NLENGTH*NWIDTH,length(hind));

%%%%%%
avglotsize=zeros(35,length(hind));
plandint=zeros(size(avglotsize));
pctdevdist=zeros(size(avglotsize));
houseperacre=zeros(size(avglotsize));
avglotsize_c=zeros(35,length(hind));
plandint_c=zeros(size(avglotsize));
pctdevdist_c=zeros(size(avglotsize));
houseperacre_c=zeros(size(avglotsize));
meandispersion=zeros(length(hind),1);
plandzone=zeros(2,length(timeset));
plandzone_c=zeros(3,length(timeset));


fullnumlots=zeros(HT,length(TSTART:TMAX),length(hind));
lotloc=cell(length(hind),1);
incomes=cell(length(hind),1);
%%%%%%
                                    
for mr=1:length(hind)
%     filename=char(fnamescell(1,hind(mr)));
%     load(eval('filename'))
%     h=strcmp(sprintf('storm_clim2_%d.mat',batchind(mr,1),...
%         batchind(mr,2)),fnamescell(1,:));
    h=strcmp(sprintf('storm_clim3_%d.mat',mr),fnamescell(1,:));
    filename=fnamescell{1,h};
    load(filename)
    
    lotloc(mr)=mat2cell(lotlocate,length(lotlocate(:,1)),2);
    ilotloc=lotloc{mr};
    
    incomes(mr)=mat2cell(cat(1,CONINFO{:,1}),length(CONINFO(:,1)),1);
    %%%%%% Density and Pct Dev w/ distance %%%%%%%%%%%%%%%%%5
%     subRlottype=Rlottype(:,:,length(testtime));
    subRlottype=reshape(LOTTYPE(:,TMAX),NLENGTH,NWIDTH);    
    subdensity=subRlottype;
    ifinddev=find(subRlottype ~= 0);
    subdensity(ifinddev)=1./z(subRlottype(ifinddev),1);
    
    startpt=0;
%     startpt=0.5195;
    endpt=startpt+0.1;
    
%     while startpt < max(max(MILEDIST))
    for ii=1:35
        % from CBD
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

        %from coast
        iring_c=find(CDIST >= startpt & CDIST < endpt);
        idev_c=(subRlottype(iring_c)~=0);
        intfarm_c=unique(setupmap(iring_c));
        ifarm_c=ismember(plandinfo(:,1),intfarm_c);
        plandint_c(ii,mr)=sum(plandinfo(ifarm_c,3).*plandinfo(ifarm_c,5))/sum(...
            plandinfo(ifarm_c,5));
        pctdevdist_c(ii,1,mr)=length(find(idev_c==1))/length(iring_c);
        houseperacre_c(ii,mr)=sum(1./z(subRlottype(iring_c(idev_c))))/length(iring_c);
%         for lt=1:HT
%         lotden(lt)=length(find(subRlottype(iring(idev))==lt))/z(lt,1);
%         end
        avglotsize_c(ii,mr)=1/(sum(1./z(subRlottype(iring_c(idev_c))))/length(iring_c(idev_c)));
        
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
    
    coastzoned=unique(setupmap(icoastal));
    coastzoned=coastzoned(coastzoned~=0);
    middlezoned=unique(setupmap(imiddle));
    middlezoned=middlezoned(middlezoned~=0);
    inlandzoned=unique(setupmap(iinland));
    inlandzoned=inlandzoned(inlandzoned~=0);
    
%     t30maps(:,:,mr)=Rlottype(:,:,length(testtime));
    t30maps(:,mr)=LOTTYPE(:,TMAX);
    subincome=zeros(NCELLS,1);
    Rbaselayer=zeros(NCELLS,1);
    for tt=1:length(testtime)
        Rbaselayer(LOTTYPE(:,testtime(tt))~=0)=1;
        subzonelots=reshape(LOTTYPE(:,testtime(tt)),NLENGTH,NWIDTH);
        subzonerents=reshape(RENT(:,testtime(tt)),NLENGTH,NWIDTH);
%         subzonelots=Rlottype(:,:,tt);
%         subzonerents=Rrent(:,:,tt);
%         ilots=find(subzonelots~=0);
        ilots=find(Rbaselayer==1);
%         ilots=cat(1,lotchoice{cat(1,lotchoice{:,4})==1,1});
        ilotloc=lotloc{mr};
        for i=1:length(ilots)
%             il=ismember(ilotloc(:,2),cat(1,lotchoice{ilots(i),2}));
            il=find(ismember(ilotloc(:,2),ilots(i))==1);
            il2=(cat(1,lotchoice{ilotloc(il,1),4})~=0);
            if isempty(find(ilotloc(il(il2),1),1))==1 | cat(1,lotchoice{ilotloc(il(il2),1),4})==0
                continue
            end
            subincome(ilotloc(il,2))=mean(cat(1,CONINFO{cat(1,lotchoice{ilotloc(il(il2),1),5}),1}));
        end
%         subincome=Rincome(:,:,tt);
        
        irents=~isnan(avgrent(:,testtime(tt)));
        rentrlts(irents,tt,mr)=avgrent(irents,testtime(tt));
        vacrlts(tt,mr)=mean(vacrate(TSTART:testtime(tt)));
%         pctutilrlts(tt,mr)=mean(pctutildiff(TSTART:testtime(tt)));
%         pctdevrlts(tt,mr)=length(find(Rbaselayer(:,:,tt)==1))/(NLENGTH*NWIDTH);
        pctdevrlts(tt,mr)=length(find(Rbaselayer==1))/NCELLS;
        pctbiglotsrlts(tt,mr)=length(find(LOTTYPE(:,testtime(tt)) > 6))/...
            length(find(Rbaselayer==1));
%         pctbiglotsrlts(tt,mr)=length(find(Rlottype(:,:,tt) > 9))/...
%             length(find(Rbaselayer(:,:,tt)==1));

        subbaselayer=reshape(Rbaselayer,NLENGTH,NWIDTH);
        unzonedpctdev(tt,mr)=length(find(subbaselayer(inotzoned)==1))/unzonedacres;
        zonedpctdev(tt,mr)=length(find(subbaselayer(izoned)==1))/zonedacres;
        coastalpctdev(tt,mr)=length(find(subbaselayer(icoastal)==1))/coastalacres;
        middlepctdev(tt,mr)=length(find(subbaselayer(imiddle)==1))/middleacres;
        inlandpctdev(tt,mr)=length(find(subbaselayer(iinland)==1))/inlandacres;
%         incomerlts(tt,mr)=consumerstats(1,testtime(tt));
        incomerlts(tt,mr)=mean(subincome(Rbaselayer==1));
        outincomerlts(tt,mr)=oldincome(testtime(tt));
%         farmerset=(1:50)';
%         idfarmnotsold=~ismember(farmerset,farmsoldinfo(:,1));
%         ifarmnotsold=ismember(setupmap,farmerset(idfarmnotsold));
        ifarmnotsold=unique(setupmap(Rbaselayer==0));
        ifarmnotsold=ifarmnotsold(ifarmnotsold~=0);
        ifarmsold=unique(setupmap(Rbaselayer==1));
        ifarmsold=ifarmsold(ifarmsold~=0);
%         subpland=Rpland(:,:,1);
        subpland=LANDINFO{3,TSTART};
        agretrlts(1,tt,mr)=mean(subpland(ismember(setupmap,ifarmsold)));
        agretrlts(2,tt,mr)=mean(subpland(ismember(setupmap,ifarmnotsold)));
%         agretrlts(1,tt,mr)=mean(Farmstats(ifarmsold,1,1));
%         agretrlts(2,tt,mr)=mean(Farmstats(ifarmnotsold,1,1));
        subcountmap=subvarmap(:,:,tt);
        ivarlayer=(Rbaselayer==1);
        subcountmap(ivarlayer)=subcountmap(ivarlayer)+1;
        subvarmap(:,:,tt)=subcountmap;

        for lt=1:HT
            ilots=find(subzonelots==lt);
            inorthlots=ismember(ilots,northindex);
            icoastlots=ismember(ilots,icoastal);
            imiddlelots=ismember(ilots,imiddle);
            if isempty(find(inorthlots,1))==0
                if tt==1
                    northrents(lt,tt,mr)=avgrent(lt,testtime(tt));
                else
                    northrents(lt,tt,mr)=mean(subzonerents(ilots(inorthlots)));
                end
            end
            if isempty(find(icoastlots,1))==0
                if tt==1
                    coastrents(lt,tt,mr)=avgrent(lt,testtime(tt));
                else
                    coastrents(lt,tt,mr)=mean(subzonerents(ilots(icoastlots)));
                end
            end
            if isempty(find(imiddlelots,1))==0
                if tt==1
                    middlerents(lt,tt,mr)=avgrent(lt,testtime(tt));
                else
                    middlerents(lt,tt,mr)=mean(subzonerents(ilots(imiddlelots)));
                end
            end
            isouthlots=ismember(ilots,southindex);
            if isempty(find(isouthlots,1))==0
                southrents(lt,tt,mr)=mean(subzonerents(ilots(isouthlots)));
            end
            iinlandlots=ismember(ilots,iinland);
            if isempty(find(iinlandlots,1))==0
                inlandrents(lt,tt,mr)=mean(subzonerents(ilots(iinlandlots)));
            end
            
            iltsize=find(subzonelots == lt);                                     
            lotsizemap(iltsize)=z(lt,1);
            
            numlotrlts(lt,tt,mr)=length(find(LOTTYPE(:,testtime(tt))==lt))/z(lt,1);
            sublotsmap=LOTTYPE(:,testtime(tt));
            unzonednumlots(lt,tt,mr)=length(find(sublotsmap(inotzoned)==lt))/z(lt,1);
            zonednumlots(lt,tt,mr)=length(find(sublotsmap(izoned)==lt))/z(lt,1);
            coastalnumlots(lt,tt,mr)=length(find(sublotsmap(icoastal)==lt))/z(lt,1);
            middlenumlots(lt,tt,mr)=length(find(sublotsmap(imiddle)==lt))/z(lt,1);
            inlandnumlots(lt,tt,mr)=length(find(sublotsmap(iinland)==lt))/z(lt,1);
            sublottypemap=lottypemap(:,:,lt,tt);
            ilottypemap=(LOTTYPE(:,testtime(tt))==lt);
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
            
            icoastincome=subincome(icoastal);
            icoastltsize=lotsizemap(icoastal);
            coastalincome(tt,mr)=mean(sum(icoastincome(icoastincome~=0).*...
                (1./icoastltsize(icoastincome~=0)))/sum(1./icoastltsize(icoastincome~=0)));
            imiddleincome=subincome(imiddle);
            imiddleltsize=lotsizemap(imiddle);
            middleincome(tt,mr)=mean(sum(imiddleincome(imiddleincome~=0).*...
                (1./imiddleltsize(imiddleincome~=0)))/sum(1./imiddleltsize(imiddleincome~=0)));
            iinlandincome=subincome(iinland);
            iinlandltsize=lotsizemap(iinland);
            inlandincome(tt,mr)=mean(sum(iinlandincome(iinlandincome~=0).*...
                (1./iinlandltsize(iinlandincome~=0)))/sum(1./iinlandltsize(iinlandincome~=0)));
        end
        
        PMAPS(:,tt,mr)=PREFMAP(:,tt);
        
    end
    totacres=length(find(ismember(setupmap,farmsoldinfo(:,1))==1));
    
    zoneavgprice(1,mr)=mean(farmsoldinfo(ismember(farmsoldinfo(:,1),outzoned),4));
    zoneavgprice(2,mr)=mean(farmsoldinfo(ismember(farmsoldinfo(:,1),inzoned),4));
    czoneavgprice(3,mr)=mean(farmsoldinfo(ismember(farmsoldinfo(:,1),coastzoned),4));
    czoneavgprice(2,mr)=mean(farmsoldinfo(ismember(farmsoldinfo(:,1),middlezoned),4));
    czoneavgprice(1,mr)=mean(farmsoldinfo(ismember(farmsoldinfo(:,1),inlandzoned),4));

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
%     zonedenrlts(:,mr)=zonedensity;
    
    Exptretrlts(:,:,mr)=Exptret/discount;
    Realretrlts(:,:,mr)=Realavgret;
    IDEALSET(:,:,mr)=idealset;
    PROFSET(:,:,mr)=profset;
    AVGBVAR(:,:,mr)=avgbrokervar;
    carcost(:,mr)=carrycost;
    
    VARLAYER(:,mr)=VARLAYER(:,mr)+BASELAYER;
    
    clear('consumerstats','vacstats','BUILDTIME','VACLAND','RENT','RETURN',...
        'LOTTYPE','BASELAYER','Rpop','Rvacrate','Rvaclots',...
        'numlt','Rleftoverpop','avgrentdynms','rentdynstats',...
        'farmsoldinfo','avgrent','avgfarmsize','stdfarmsize','DELTA',...
        'survivalrate','LOCWGHT','REGWGHT','PCTSEARCH','zeta','HIBETA',...
        'MIDBETA','LOWBETA','POPGROW','ccost','newhouseset','bidtot',...
        'meandisp','maxdevdist',...
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
% 
% %%%%%% Build VARMAP for each time mark %%%%%
% agretrlts(1,1,:)=0;
% for tt=1:length(testtime)
%     VARMAP(:,:,tt)=subvarmap(:,:,tt)./length(hind);
% 
% 
% 
%     %%%%%% Load last run's file %%%%%%%%%
%     % load resultsLOC30_160210_20.mat
% 
%     % VARLAYERSAVE(:,:)=VARLAYER;
%     % varmap=VARLAYERSAVE/max(max(VARLAYER));
%     ivarthresh=(VARMAP(:,:,tt) >= avgthresh(tt));
%     locmat=find(ivarthresh==1);
% 
%     lottypemap(:,:,:,tt)=lottypemap(:,:,:,tt)./length(hind);
%     [maxprob,imaxprob]=max(lottypemap(:,:,:,tt),[],3);
% 
% 
%     for lt=1:HT
%         singleltmap=lottypemap(:,:,lt,tt);
%         iltzero=(imaxprob==lt);
%         singleltmap(iltzero)=0;
%         secprobmap(:,:,lt,tt)=singleltmap;
%     end
%     [secprob,isecprob]=max(secprobmap(:,:,:,tt),[],3);
% 
%     confprob=maxprob-secprob;
% %     PRODMAP=Rpland(:,:,1);
%     %%%%% Test representativeness of varmap
%     testpctdev(tt)=length(find(ivarthresh==1))/(NLENGTH*NWIDTH);
%     imapcover=find(ivarthresh==1);
%     testmap=zeros(NLENGTH,NWIDTH);
%     testmap(imapcover)=1;
%     altprop=zeros(length(find(imapcover)),1);
%     for ic=1:length(imapcover)
%         [row,col]=ind2sub([NLENGTH NWIDTH],imapcover(ic));
%         updir=max(row-2,1);
%         dndir=min(row+2,NLENGTH);
%         lfdir=max(col-2,1);
%         rtdir=min(col+2,NWIDTH);
% 
% 
%         altprop(ic)=length(find(testmap(updir:dndir,lfdir:rtdir)==0))/...
%             (length(updir:dndir)*length(lfdir:rtdir));
%     end
%     subAVGMAP=zeros(NLENGTH,NWIDTH);
%     if tt > 1
%         idevthere=(AVGMAP(:,:,tt-1)~=0);
%         indthresh=find(ivarthresh==1);
%         inddev=find(idevthere==1);
%         idevadd=~ismember(indthresh,inddev);
%         subAVGMAP=AVGMAP(:,:,tt-1);
%         subAVGMAP(indthresh(idevadd))=imaxprob(indthresh(idevadd));
%         AVGMAP(:,:,tt)=subAVGMAP;
%     else
%         subAVGMAP(ivarthresh)=imaxprob(ivarthresh);
%         AVGMAP(:,:,tt)=subAVGMAP;
%     end
%     testmeandisp=sum(altprop)/length(imapcover);
%     subCONFMAP=zeros(NLENGTH,NWIDTH);
%     subCONFMAP(ivarthresh)=confprob(ivarthresh);
%     CONFMAP(:,:,tt)=subCONFMAP;
%     
%     timePMAPS=zeros(NLENGTH*NWIDTH,length(hind));
%     timePMAPS(:,:)=PMAPS(:,tt,:);
%     for il=1:length(locmat)
% %         [irow,icol]=ind2sub([NLENGTH NWIDTH],locmat(il));
% %         isnotzero=(timePMAPS(irow,icol,:)~=0);
%         isnotzero=(timePMAPS(locmat(il),:)~=0);
%         avgPREFMAP(locmat(il),tt)=mean(timePMAPS(locmat(il),isnotzero));
%     end
%     ifindnan=find(isnan(avgPREFMAP(:,tt))==1);  
% end

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
avgidealset=mean(IDEALSET(:,testtime,:),3);
diffset=mean(IDEALSET(:,testtime,:)-numlotrlts,3);

avgprofset=mean(PROFSET(:,testtime,:),3);
diffprofset=mean(PROFSET(:,testtime,:)-numlotrlts,3);

avgvar=mean(AVGBVAR(:,10:30,:),3);

%%
%     
% figure(1)
% surf(AVGMAP(:,:,tt))
% axis ij
% view(0,90)
% title('Average Landscape Outcome, t=30')
% set(gca,'clim',[0 HT])
% colorbar
% 
% figure(2)
% surf(VARMAP(:,:,5))
% axis ij
% view(0,90)
% title('Probability of Development, t=30')
% colorbar
% 
% figure(3)
% pcolor(avgPREFMAP(:,tt))
% axis ij
% title('Average Consumer Preference for Proximity to Coast')
% set(gca,'clim',[0 0.3])
% colorbar
% 
% % % 
% % % figure(2)
% % surf(CONFMAP(:,:,2))
% % axis ij
% % view(0,90)
% % title('Likelihood of Average Landscape Outcome, No Zoning, t=15')
% % colorbar


%% Generate Statistics %%

subavglotsize=zeros(size(avglotsize));
subavglotsize(~isnan(avglotsize))=avglotsize(~isnan(avglotsize));
subhouseperacre=zeros(size(houseperacre));
subhouseperacre(~isnan(houseperacre))=houseperacre(~isnan(houseperacre));
subplandint=zeros(size(avglotsize));
subplandint(~isnan(plandint))=plandint(~isnan(plandint));

subavglotsize_c=zeros(size(avglotsize_c));
subavglotsize_c(~isnan(avglotsize_c))=avglotsize_c(~isnan(avglotsize_c));
subhouseperacre_c=zeros(size(houseperacre_c));
subhouseperacre_c(~isnan(houseperacre_c))=houseperacre_c(~isnan(houseperacre_c));
subplandint_c=zeros(size(avglotsize_c));
subplandint_c(~isnan(plandint_c))=plandint_c(~isnan(plandint_c));

northfarmers=unique(setupmap(ZONEMAP==0));
northfarmers=northfarmers(northfarmers~=0);
southfarmers=unique(setupmap(ZONEMAP==1));
inorthfarmers=ismember(plandinfo(:,1),northfarmers);
isouthfarmers=ismember(plandinfo(:,1),southfarmers);

coastfarmers=unique(setupmap(CZONEMAP==3));
coastfarmers=coastfarmers(coastfarmers~=0);
middlefarmers=unique(setupmap(CZONEMAP==2));
middlefarmers=middlefarmers(middlefarmers~=0);
inlandfarmers=unique(setupmap(CZONEMAP==1));
icoastfarmers=ismember(plandinfo(:,1),coastfarmers);
imiddlefarmers=ismember(plandinfo(:,1),middlefarmers);
iinlandfarmers=ismember(plandinfo(:,1),inlandfarmers);


for imr=1:35
    inonzerolot=(subavglotsize(imr,:) ~= 0);
    inonzerohouse=(subhouseperacre(imr,:) ~= 0);
    inonzeroplandint=(subplandint(imr,:) ~= 0);
    meanlotsize(imr)=mean(subavglotsize(imr,inonzerolot));
    meanhouseperacre(imr)=mean(subhouseperacre(imr,inonzerohouse));
    meanplandint(imr)=mean(subplandint(imr,inonzeroplandint));
    
    inonzerolot_c=(subavglotsize_c(imr,:) ~= 0);
    inonzerohouse_c=(subhouseperacre_c(imr,:) ~= 0);
    inonzeroplandint_c=(subplandint(imr,:) ~= 0);
    meanlotsize_c(imr)=mean(subavglotsize_c(imr,inonzerolot_c));
    meanhouseperacre_c(imr)=mean(subhouseperacre_c(imr,inonzerohouse_c));
    meanplandint_c(imr)=mean(subplandint_c(imr,inonzeroplandint_c));
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
   
   icfarm=icoastfarmers.*itimeset;
   imidfarm=imiddlefarmers.*itimeset;
   iinfarm=iinlandfarmers.*itimeset;
   plandzone_c(3,pt)=mean(plandinfo(icfarm==1,3));
   plandzone_c(2,pt)=mean(plandinfo(imidfarm==1,3));
   plandzone_c(1,pt)=mean(plandinfo(iinfarm==1,3));
   
%    avgplandtime(pt,2)=mean(plandinfo(itimeset,3).*plandinfo(itimeset,5));
   avgplandtime(pt,2)=sum(plandinfo(itimeset,3).*plandinfo(itimeset,5))/...
       sum(plandinfo(itimeset,5));
   
   for lt=1:HT
       subrets(1:length(hind))=retrlts(lt,pt,:);
       iposrets=(subrets~=0);
%        subhtvac(1:length(hind))=htutilvacrlts(lt,pt,:);
%        iposhtv=(subhtvac~=0);
%        subhtnovac(1:length(hind))=htutilnovacrlts(lt,pt,:);
%        iposhtnv=(subhtnovac~=0);
       subhtinc(1:length(hind))=htincomerlts(lt,pt,:);
       iposhti=(subhtinc~=0);
       subbids(1:length(hind))=bidsharerlts(lt,pt,:);
       iposbids=(subbids~=0);
       
       subretlots(1:length(hind))=fullnumlots(lt,pt+1,:);
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
plandzone_c(isnan(plandzone_c))=0;
% for pd=1:50
for pd=1:length(planddist)
   avgplanddist(pd,1)=planddist(pd).*0.0395;
   avgplanddist(pd,2)=mean(plandinfo(plandinfo(:,4)==planddist(pd),3));
end

carrycostrlts(1,:)=mean(carcost,2)';
for tt=1:length(testtime)
    for lt=1:HT
        subnorthrents(1:length(hind))=northrents(lt,tt,:);
        iposnorthrents=(subnorthrents~=0);
        subsouthrents(1:length(hind))=southrents(lt,tt,:);
        ipossouthrents=(subsouthrents~=0);
        subnorthlots(1:length(hind))=unzonednumlots(lt,tt,:);
        iposnorthlots=(subnorthlots~=0);
        subsouthlots(1:length(hind))=zonednumlots(lt,tt,:);
        ipossouthlots=(subsouthlots~=0);
        
        subcoastrents(1:length(hind))=coastrents(lt,tt,:);
        iposcoastrents=(subcoastrents~=0);
        submiddlerents(1:length(hind))=middlerents(lt,tt,:);
        iposmiddlerents=(submiddlerents~=0);
        subinlandrents(1:length(hind))=inlandrents(lt,tt,:);
        iposinlandrents=(subinlandrents~=0);
        subcoastlots(1:length(hind))=coastalnumlots(lt,tt,:);
        iposcoastlots=(subcoastlots~=0);
        submiddlelots(1:length(hind))=middlenumlots(lt,tt,:);
        iposmiddlelots=(submiddlelots~=0);
        subinlandlots(1:length(hind))=inlandnumlots(lt,tt,:);
        iposinlandlots=(subinlandlots~=0);
        
        subrents(1:length(hind))=rentrlts(lt,tt,:);
        iposrents=(subrents~=0);
        sublots(1:length(hind))=numlotrlts(lt,tt,:);
        iposlots=(sublots~=0);
        
        [nrentmu(lt,tt),nrentsigma(lt,tt),nrentmuci(:,lt,tt),nrentsigmaci(:,lt,tt)]=normfit(subnorthrents(iposnorthrents));
        nrentmu(lt,tt)=sum(subnorthrents(iposnorthrents).*subnorthlots(iposnorthrents))/sum(subnorthlots(iposnorthrents));
        nrentse(lt,tt)=nrentsigma(lt,tt)/sqrt(length(subnorthrents(iposnorthrents)));
        
        [srentmu(lt,tt),srentsigma(lt,tt),srentmuci(:,lt,tt),srentsigmaci(:,lt,tt)]=normfit(subsouthrents(ipossouthrents));
        srentmu(lt,tt)=sum(subsouthrents(ipossouthrents).*subsouthlots(ipossouthrents))/sum(subsouthlots(ipossouthrents));
        srentse(lt,tt)=srentsigma(lt,tt)/sqrt(length(subsouthrents(ipossouthrents)));
        
        [crentmu(lt,tt),crentsigma(lt,tt),crentmuci(:,lt,tt),crentsigmaci(:,lt,tt)]=normfit(subcoastrents(iposcoastrents));
        crentmu(lt,tt)=sum(subcoastrents(iposcoastrents).*subcoastlots(iposcoastrents))/sum(subcoastlots(iposcoastrents));
        crentse(lt,tt)=crentsigma(lt,tt)/sqrt(length(subcoastrents(iposcoastrents)));
        
        [midrentmu(lt,tt),midrentsigma(lt,tt),midrentmuci(:,lt,tt),midrentsigmaci(:,lt,tt)]=normfit(submiddlerents(iposmiddlerents));
        midrentmu(lt,tt)=sum(submiddlerents(iposmiddlerents).*submiddlelots(iposmiddlerents))/sum(submiddlelots(iposmiddlerents));
        midrentse(lt,tt)=midrentsigma(lt,tt)/sqrt(length(submiddlerents(iposmiddlerents)));
        
        [inrentmu(lt,tt),inrentsigma(lt,tt),inrentmuci(:,lt,tt),inrentsigmaci(:,lt,tt)]=normfit(subinlandrents(iposinlandrents));
        inrentmu(lt,tt)=sum(subinlandrents(iposinlandrents).*subinlandlots(iposinlandrents))/sum(subinlandlots(iposinlandrents));
        inrentse(lt,tt)=inrentsigma(lt,tt)/sqrt(length(subinlandrents(iposinlandrents)));
        
        if isempty(find(iposrents,1))==1
            continue
        end
        [rentmu(lt,tt),rentsigma(lt,tt),rentmuci(:,lt,tt),rentsigmaci(:,lt,tt)]=normfit(subrents(iposrents));
        rentmu(lt,tt)=sum(subrents(iposrents).*sublots(iposrents))/sum(sublots(iposrents));
        rentse(lt,tt)=rentsigma(lt,tt)/sqrt(length(subrents(iposrents)));
%         [lotmu(lt,tt),lotsigma(lt,tt),lotmuci(:,lt,tt),lotsigmaci(:,lt,tt)]=normfit(sublots(iposlots));
        [lotmu(lt,tt),lotsigma(lt,tt),lotmuci(:,lt,tt),lotsigmaci(:,lt,tt)]=normfit(sublots);        
        lotse(lt,tt)=lotsigma(lt,tt)/sqrt(length(sublots));
        subunzonednumlots(1:length(hind))=unzonednumlots(lt,tt,:);
        subzonednumlots(1:length(hind))=zonednumlots(lt,tt,:);
        [unzonedlotmu(lt,tt),unzonedlotsigma(lt,tt),unzonedlotmuci(:,lt,tt),unzonedlotsigmaci(:,lt,tt)]=normfit(subunzonednumlots);        
        unzonedlotse(lt,tt)=unzonedlotsigma(lt,tt)/sqrt(length(subunzonednumlots));
        [zonedlotmu(lt,tt),zonedlotsigma(lt,tt),zonedlotmuci(:,lt,tt),zonedlotsigmaci(:,lt,tt)]=normfit(subzonednumlots);        
        zonedlotse(lt,tt)=zonedlotsigma(lt,tt)/sqrt(length(subzonednumlots));
        
        subcoastalnumlots(1:length(hind))=coastalnumlots(lt,tt,:);
        submiddlenumlots(1:length(hind))=middlenumlots(lt,tt,:);
        subinlandnumlots(1:length(hind))=inlandnumlots(lt,tt,:);
        [coastlotmu(lt,tt),coastlotsigma(lt,tt),coastlotmuci(:,lt,tt),coastlotsigmaci(:,lt,tt)]=normfit(subcoastalnumlots);        
        coastlotse(lt,tt)=coastlotsigma(lt,tt)/sqrt(length(coastalnumlots));
        [middlelotmu(lt,tt),middlelotsigma(lt,tt),middlelotmuci(:,lt,tt),middlelotsigmaci(:,lt,tt)]=normfit(submiddlenumlots);        
        middlelotse(lt,tt)=middlelotsigma(lt,tt)/sqrt(length(middlenumlots));
        [inlandlotmu(lt,tt),inlandlotsigma(lt,tt),inlandlotmuci(:,lt,tt),inlandlotsigmaci(:,lt,tt)]=normfit(subinlandnumlots);        
        inlandlotse(lt,tt)=inlandlotsigma(lt,tt)/sqrt(length(inlandnumlots));
    end
    for zl=1:length(zonedenrlts(:,1))
        iposzone=~isnan(zonedenrlts(zl,:));
        [zonedenmu(zl),zonedensigma(zl),zonedenmuci(:,zl),zonedensigmaci(:,zl)]=normfit(zonedenrlts(zl,iposzone));
        zonedense(zl)=zonedensigma(zl)/sqrt(length(hind));
    end
    for zl=1:length(czonedenrlts(:,1))    
        iposczone=~isnan(czonedenrlts(zl,:));
        [czonedenmu(zl),czonedensigma(zl),czonedenmuci(:,zl),czonedensigmaci(:,zl)]=normfit(czonedenrlts(zl,iposczone));
        czonedense(zl)=czonedensigma(zl)/sqrt(length(hind));
    end
    

    % keyboard    %may need to transpose, check size of rentmu
    nrentstats(:,:,tt)=max([nrentmu(:,tt),nrentsigma(:,tt),nrentmuci(1,:,tt)',nrentmuci(2,:,tt)',nrentse(:,tt)],0);
    srentstats(:,:,tt)=max([srentmu(:,tt),srentsigma(:,tt),srentmuci(1,:,tt)',srentmuci(2,:,tt)',srentse(:,tt)],0);
       
    crentstats(:,:,tt)=max([crentmu(:,tt),crentsigma(:,tt),crentmuci(1,:,tt)',crentmuci(2,:,tt)',crentse(:,tt)],0);
    midrentstats(:,:,tt)=max([midrentmu(:,tt),midrentsigma(:,tt),midrentmuci(1,:,tt)',midrentmuci(2,:,tt)',midrentse(:,tt)],0);
    inrentstats(:,:,tt)=max([inrentmu(:,tt),inrentsigma(:,tt),inrentmuci(1,:,tt)',inrentmuci(2,:,tt)',inrentse(:,tt)],0);

    
    rentstats(:,:,tt)=max([rentmu(:,tt),rentsigma(:,tt),rentmuci(1,:,tt)',rentmuci(2,:,tt)',rentse(:,tt)],0);
    lotstats(:,:,tt)=max([lotmu(:,tt),lotsigma(:,tt),lotmuci(1,:,tt)',lotmuci(2,:,tt)',lotse(:,tt)],0);
    
    
    unzonedlotstats(:,:,tt)=max([unzonedlotmu(:,tt),unzonedlotsigma(:,tt),unzonedlotmuci(1,:,tt)',unzonedlotmuci(2,:,tt)',unzonedlotse(:,tt)],0);
    zonedlotstats(:,:,tt)=max([zonedlotmu(:,tt),zonedlotsigma(:,tt),zonedlotmuci(1,:,tt)',zonedlotmuci(2,:,tt)',zonedlotse(:,tt)],0);
    
    coastlotstats(:,:,tt)=max([coastlotmu(:,tt),coastlotsigma(:,tt),coastlotmuci(1,:,tt)',coastlotmuci(2,:,tt)',coastlotse(:,tt)],0);
    middlelotstats(:,:,tt)=max([middlelotmu(:,tt),middlelotsigma(:,tt),middlelotmuci(1,:,tt)',middlelotmuci(2,:,tt)',middlelotse(:,tt)],0);
    inlandlotstats(:,:,tt)=max([inlandlotmu(:,tt),inlandlotsigma(:,tt),inlandlotmuci(1,:,tt)',inlandlotmuci(2,:,tt)',inlandlotse(:,tt)],0);

    [lotsummu(tt),lotsumsigma(tt),lotsummuci(:,tt),lotsumsigmaci(:,tt)]=normfit(lotsumrlts(tt,:));
    lotsumse(tt)=lotsumsigma(tt)/sqrt(length(hind));
    lotsumstats(tt,:)=max([lotsummu(tt),lotsumsigma(tt),lotsummuci(1,tt),lotsummuci(2,tt),lotsumse(tt)],0);

    [vacmu(tt),vacsigma(tt),vacmuci(:,tt),vacsigmaci(:,tt)]=normfit(vacrlts(tt,:));
    vacse(tt)=vacsigma(tt)/sqrt(length(hind));
    vacstats(tt,:)=max([vacmu(tt),vacsigma(tt),vacmuci(1,tt),vacmuci(2,tt),vacse(tt)],0);
    
    [incomemu(tt),incomesigma(tt),incomemuci(:,tt),incomesigmaci(:,tt)]=normfit(incomerlts(tt,:));
    incomese(tt)=incomesigma(tt)/sqrt(length(hind));
    incomestats(tt,:)=max([incomemu(tt),incomesigma(tt),incomemuci(1,tt),incomemuci(2,tt),incomese(tt)],0);
    
    [northincomemu(tt),northincomesigma(tt),northincomemuci(:,tt),northincomesigmaci(:,tt)]=normfit(northincome(tt,:));
    northincomese(tt)=northincomesigma(tt)/sqrt(length(hind));
    northincomestats(tt,:)=max([northincomemu(tt),northincomesigma(tt),northincomemuci(1,tt),northincomemuci(2,tt),northincomese(tt)],0);
    
    [southincomemu(tt),southincomesigma(tt),southincomemuci(:,tt),southincomesigmaci(:,tt)]=normfit(southincome(tt,:));
    southincomese(tt)=southincomesigma(tt)/sqrt(length(hind));
    southincomestats(tt,:)=max([southincomemu(tt),southincomesigma(tt),southincomemuci(1,tt),southincomemuci(2,tt),southincomese(tt)],0);
    
    [coastincomemu(tt),coastincomesigma(tt),coastincomemuci(:,tt),coastincomesigmaci(:,tt)]=normfit(coastalincome(tt,:));
    coastincomese(tt)=coastincomesigma(tt)/sqrt(length(hind));
    coastincomestats(tt,:)=max([coastincomemu(tt),coastincomesigma(tt),coastincomemuci(1,tt),coastincomemuci(2,tt),coastincomese(tt)],0);
    
    [middleincomemu(tt),middleincomesigma(tt),middleincomemuci(:,tt),middleincomesigmaci(:,tt)]=normfit(middleincome(tt,:));
    middleincomese(tt)=middleincomesigma(tt)/sqrt(length(hind));
    middleincomestats(tt,:)=max([middleincomemu(tt),middleincomesigma(tt),middleincomemuci(1,tt),middleincomemuci(2,tt),middleincomese(tt)],0);
    
    [inlandincomemu(tt),inlandincomesigma(tt),inlandincomemuci(:,tt),inlandincomesigmaci(:,tt)]=normfit(inlandincome(tt,:));
    inlandincomese(tt)=inlandincomesigma(tt)/sqrt(length(hind));
    inlandincomestats(tt,:)=max([inlandincomemu(tt),inlandincomesigma(tt),inlandincomemuci(1,tt),inlandincomemuci(2,tt),inlandincomese(tt)],0);

    [outincomemu(tt),outincomesigma(tt),outincomemuci(:,tt),outincomesigmaci(:,tt)]=normfit(outincomerlts(tt,:));
    outincomese(tt)=outincomesigma(tt)/sqrt(length(hind));
    outincomestats(tt,:)=max([outincomemu(tt),outincomesigma(tt),outincomemuci(1,tt),outincomemuci(2,tt),outincomese(tt)],0);
    
    [zonedpctdevmu(tt),zonedpctdevsigma(tt),zonedpctdevmuci(:,tt),zonedpctdevsigmaci(:,tt)]=normfit(zonedpctdev(tt,:));
    zonedpctdevse(tt)=zonedpctdevsigma(tt)/sqrt(length(hind));
    zonedpctdevstats(tt,:)=max([zonedpctdevmu(tt),zonedpctdevsigma(tt),zonedpctdevmuci(1,tt),zonedpctdevmuci(2,tt),zonedpctdevse(tt)],0);
    
    [unzonedpctdevmu(tt),unzonedpctdevsigma(tt),unzonedpctdevmuci(:,tt),unzonedpctdevsigmaci(:,tt)]=normfit(unzonedpctdev(tt,:));
    unzonedpctdevse(tt)=unzonedpctdevsigma(tt)/sqrt(length(hind));
    unzonedpctdevstats(tt,:)=max([unzonedpctdevmu(tt),unzonedpctdevsigma(tt),unzonedpctdevmuci(1,tt),unzonedpctdevmuci(2,tt),unzonedpctdevse(tt)],0);
    
    [coastpctdevmu(tt),coastpctdevsigma(tt),coastpctdevmuci(:,tt),coastpctdevsigmaci(:,tt)]=normfit(coastalpctdev(tt,:));
    coastpctdevse(tt)=coastpctdevsigma(tt)/sqrt(length(hind));
    coastpctdevstats(tt,:)=max([coastpctdevmu(tt),coastpctdevsigma(tt),coastpctdevmuci(1,tt),coastpctdevmuci(2,tt),coastpctdevse(tt)],0);
    
    [middlepctdevmu(tt),middlepctdevsigma(tt),middlepctdevmuci(:,tt),middlepctdevsigmaci(:,tt)]=normfit(middlepctdev(tt,:));
    middlepctdevse(tt)=middlepctdevsigma(tt)/sqrt(length(hind));
    middlepctdevstats(tt,:)=max([middlepctdevmu(tt),middlepctdevsigma(tt),middlepctdevmuci(1,tt),middlepctdevmuci(2,tt),middlepctdevse(tt)],0);
   
    [inlandpctdevmu(tt),inlandpctdevsigma(tt),inlandpctdevmuci(:,tt),inlandpctdevsigmaci(:,tt)]=normfit(inlandpctdev(tt,:));
    inlandpctdevse(tt)=inlandpctdevsigma(tt)/sqrt(length(hind));
    inlandpctdevstats(tt,:)=max([inlandpctdevmu(tt),inlandpctdevsigma(tt),inlandpctdevmuci(1,tt),inlandpctdevmuci(2,tt),inlandpctdevse(tt)],0);

    [utilmu(tt),utilsigma(tt),utilmuci(:,tt),utilsigmaci(:,tt)]=normfit(pctutilrlts(tt,:));
    utilse(tt)=utilsigma(tt)/sqrt(length(hind));
    utilstats(tt,:)=max([utilmu(tt),utilsigma(tt),utilmuci(1,tt),utilmuci(2,tt),utilse(tt)],0);
    
    [devmu(tt),devsigma(tt),devmuci(:,tt),devsigmaci(:,tt)]=normfit(pctdevrlts(tt,:));
    devse(tt)=devsigma(tt)/sqrt(length(hind));
    devstats(tt,:)=max([devmu(tt),devsigma(tt),devmuci(1,tt),devmuci(2,tt),devse(tt)],0);
    
    [biglotsmu(tt),biglotssigma(tt),biglotsmuci(:,tt),biglotssigmaci(:,tt)]=normfit(pctbiglotsrlts(tt,:));
    biglotsse(tt)=biglotssigma(tt)/sqrt(length(hind));
    biglotsstats(tt,:)=max([biglotsmu(tt),biglotssigma(tt),biglotsmuci(1,tt),biglotsmuci(2,tt),biglotsse(tt)],0);
    
    subagretsold(1:length(hind))=agretrlts(1,tt,:);
    [agretsoldmu(tt),agretsoldsigma(tt),agretsoldmuci(:,tt),agretsoldsigmaci(:,tt)]=normfit(subagretsold);
    agretsoldse(tt)=agretsoldsigma(tt)/sqrt(length(hind));
    agretsoldstats(tt,:)=max([agretsoldmu(tt),agretsoldsigma(tt),agretsoldmuci(1,tt),agretsoldmuci(2,tt),agretsoldse(tt)],0);
    
    subagretfarm(1:length(hind))=agretrlts(2,tt,:);
    [agretfarmmu(tt),agretfarmsigma(tt),agretfarmmuci(:,tt),agretfarmsigmaci(:,tt)]=normfit(subagretfarm);
    agretfarmse(tt)=agretfarmsigma(tt)/sqrt(length(hind));
    agretfarmstats(tt,:)=max([agretfarmmu(tt),agretfarmsigma(tt),agretfarmmuci(1,tt),agretfarmmuci(2,tt),agretfarmse(tt)],0);
    
end

zonedenstats=max([zonedenmu,zonedensigma,zonedenmuci(1,:)',zonedenmuci(2,:)',zonedense],0);
[dispmu,dispsigma,dispmuci,dispsigmaci]=normfit(disprlts);
dispse=dispsigma/sqrt(length(hind));
dispstats=max([dispmu,dispsigma,dispmuci(1),dispmuci(2),dispse],0);
[distmu,distsigma,distmuci,distsigmaci]=normfit(maxdistrlts);
distse=distsigma/sqrt(length(hind));
diststats=max([distmu,distsigma,distmuci(1),distmuci(2),distse],0);

[coastpricemu,coastpricesigma,coastpricemuci,coastpricesigmaci]=normfit(czoneavgprice(3,:));
coastpricese=coastpricesigma/sqrt(length(hind));
czoneavgpricestats(3,:)=max([coastpricemu,coastpricesigma,coastpricemuci(1),coastpricemuci(2),coastpricese],0);

[middlepricemu,middlepricesigma,middlepricemuci,middlepricesigmaci]=normfit(czoneavgprice(2,:));
middlepricese=middlepricesigma/sqrt(length(hind));
czoneavgpricestats(2,:)=max([middlepricemu,middlepricesigma,middlepricemuci(1),middlepricemuci(2),middlepricese],0);

[inlandpricemu,inlandpricesigma,inlandpricemuci,inlandpricesigmaci]=normfit(czoneavgprice(1,:));
inlandpricese=inlandpricesigma/sqrt(length(hind));
czoneavgpricestats(1,:)=max([inlandpricemu,inlandpricesigma,inlandpricemuci(1),inlandpricemuci(2),inlandpricese],0);
    
[unzonedpricemu,unzonedpricesigma,unzonedpricemuci,unzonedpricesigmaci]=normfit(zoneavgprice(1,:));
unzonedpricese=unzonedpricesigma/sqrt(length(hind));
zoneavgpricestats(1,:)=max([unzonedpricemu,unzonedpricesigma,unzonedpricemuci(1),unzonedpricemuci(2),unzonedpricese],0);
[zonedpricemu,zonedpricesigma,zonedpricemuci,zonedpricesigmaci]=normfit(zoneavgprice(2,:));
zonedpricese=zonedpricesigma/sqrt(length(hind));
zoneavgpricestats(2,:)=max([zonedpricemu,zonedpricesigma,zonedpricemuci(1),zonedpricemuci(2),zonedpricese],0);

[overplandmu,overplandsigma,overplandmuci,overplandsigmaci]=normfit(overpland(1,:));
overplandse=overplandsigma/sqrt(length(hind));
overplandstats(1,:)=max([overplandmu,overplandsigma,overplandmuci(1),overplandmuci(2),overplandse],0);
%     [psoldmu,psoldsigma,psoldmuci,psoldsigmaci]=normfit(plandrlts(1,:));
%     psoldse=psoldsigma/sqrt(length(hind));
%     psoldstats=max([psoldmu,psoldsigma,psoldmuci(1),psoldmuci(2),psoldse],0);
%     [pfarmmu,pfarmsigma,pfarmmuci,pfarmsigmaci]=normfit(plandrlts(2,:));
%     pfarmse=pfarmsigma/sqrt(length(hind));
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
distlabel_c={'Distance from Coast'};
incomelabel={'Mean Income of Occupants'};
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
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\results\alt_storm_climate
resultsfile=('Results_CHALMS_alt_storm_clim3.xlsx');
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
xlswrite(resultsfile,lotlabel,'Time Step Stats','A12');
xlswrite(resultsfile,lottypeset,'Time Step Stats','B3');
xlswrite(resultsfile,lottypeset,'Time Step Stats','B12');
xlswrite(resultsfile,lotsumlabel,'Time Step Stats','B20');
xlswrite(resultsfile,statlabel,'Time Step Stats','C22');
xlswrite(resultsfile,distlabel,'Time Step Stats','B22');
xlswrite(resultsfile,zonedenlabel,'Time Step Stats','A22');
xlswrite(resultsfile,plandtimelabel,'Time Step Stats','B31');
xlswrite(resultsfile,tlabel,'Time Step Stats','A32');
xlswrite(resultsfile,pllabel,'Time Step Stats','A33');
xlswrite(resultsfile,planddistlabel,'Time Step Stats','B35');
xlswrite(resultsfile,dlabel,'Time Step Stats','A36');
xlswrite(resultsfile,pllabel,'Time Step Stats','A37');

xlswrite(resultsfile,rentstats(:,:,1),'Time Step Stats','C3');
xlswrite(resultsfile,rentstats(:,:,2),'Time Step Stats','I3');
xlswrite(resultsfile,rentstats(:,:,3),'Time Step Stats','O3');
xlswrite(resultsfile,rentstats(:,:,4),'Time Step Stats','U3');
xlswrite(resultsfile,rentstats(:,:,5),'Time Step Stats','AA3');
xlswrite(resultsfile,lotstats(:,:,1),'Time Step Stats','C12');
xlswrite(resultsfile,lotstats(:,:,2),'Time Step Stats','I12');
xlswrite(resultsfile,lotstats(:,:,3),'Time Step Stats','O12');
xlswrite(resultsfile,lotstats(:,:,4),'Time Step Stats','U12');
xlswrite(resultsfile,lotstats(:,:,5),'Time Step Stats','AA12');
xlswrite(resultsfile,lotsumstats(1,:),'Time Step Stats','C20');
xlswrite(resultsfile,lotsumstats(2,:),'Time Step Stats','I20');
xlswrite(resultsfile,lotsumstats(3,:),'Time Step Stats','O20');
xlswrite(resultsfile,lotsumstats(4,:),'Time Step Stats','U20');
xlswrite(resultsfile,lotsumstats(5,:),'Time Step Stats','AA20');

% xlswrite(resultsfile,zonedenstats,'Time Step Stats','C33');
% xlswrite(resultsfile,zonedist,'Time Step Stats','B33');
xlswrite(resultsfile,maxdistlabel,'Time Step Stats','A22');
xlswrite(resultsfile,statlabel,'Time Step Stats','A23');
xlswrite(resultsfile,diststats,'Time Step Stats','A24');
xlswrite(resultsfile,displabel,'Time Step Stats','A26');
xlswrite(resultsfile,statlabel,'Time Step Stats','A27');
xlswrite(resultsfile,dispstats,'Time Step Stats','A28');

xlswrite(resultsfile,avgplandtime','Time Step Stats','B32');
xlswrite(resultsfile,avgplanddist','Time Step Stats','B36');
xlswrite(resultsfile,carrycostlabel,'Time Step Stats','A34');
xlswrite(resultsfile,carrycostrlts(1,TSTART+1:TMAX),'Time Step Stats','B34');

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

xlswrite(resultsfile,testtime','Landscape Results','V3');
xlswrite(resultsfile,coastalpctdevlabel,'Landscape Results','W1');
xlswrite(resultsfile,statlabel,'Landscape Results','W2');
xlswrite(resultsfile,coastpctdevstats,'Landscape Results','W3');
xlswrite(resultsfile,testtime','Landscape Results','V11');
xlswrite(resultsfile,middlepctdevlabel,'Landscape Results','W9');
xlswrite(resultsfile,statlabel,'Landscape Results','W10');
xlswrite(resultsfile,middlepctdevstats,'Landscape Results','W11');
xlswrite(resultsfile,testtime','Landscape Results','V19');
xlswrite(resultsfile,inlandpctdevlabel,'Landscape Results','W17');
xlswrite(resultsfile,statlabel,'Landscape Results','W18');
xlswrite(resultsfile,inlandpctdevstats,'Landscape Results','W19');

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

xlswrite(resultsfile,coastpricelabel,'Landscape Results','B29');
xlswrite(resultsfile,statlabel,'Landscape Results','B30');
xlswrite(resultsfile,inlandlabel,'Landscape Results','A31');
xlswrite(resultsfile,middlelabel,'Landscape Results','A32');
xlswrite(resultsfile,coastlabel,'Landscape Results','A33');
xlswrite(resultsfile,czoneavgpricestats,'Landscape Results','B31');

xlswrite(resultsfile,overplandlabel,'Landscape Results','B34');
xlswrite(resultsfile,statlabel,'Landscape Results','B35');
xlswrite(resultsfile,overplandstats,'Landscape Results','B36');

xlswrite(resultsfile,northincomelabel,'Landscape Results','P1');
xlswrite(resultsfile,statlabel,'Landscape Results','P2');
xlswrite(resultsfile,testtime','Landscape Results','O3');
xlswrite(resultsfile,northincomestats,'Landscape Results','P3');
xlswrite(resultsfile,southincomelabel,'Landscape Results','P9');
xlswrite(resultsfile,statlabel,'Landscape Results','P10');
xlswrite(resultsfile,testtime','Landscape Results','O11');
xlswrite(resultsfile,southincomestats,'Landscape Results','P11');

xlswrite(resultsfile,coastincomelabel,'Landscape Results','P18');
xlswrite(resultsfile,statlabel,'Landscape Results','P19');
xlswrite(resultsfile,testtime','Landscape Results','O20');
xlswrite(resultsfile,coastincomestats,'Landscape Results','P20');
xlswrite(resultsfile,middleincomelabel,'Landscape Results','P26');
xlswrite(resultsfile,statlabel,'Landscape Results','P27');
xlswrite(resultsfile,testtime','Landscape Results','O28');
xlswrite(resultsfile,middleincomestats,'Landscape Results','P28');
xlswrite(resultsfile,inlandincomelabel,'Landscape Results','P34');
xlswrite(resultsfile,statlabel,'Landscape Results','P35');
xlswrite(resultsfile,testtime','Landscape Results','O36');
xlswrite(resultsfile,inlandincomestats,'Landscape Results','P36');

xlswrite(resultsfile,distlabel,'Landscape Results','B41');
xlswrite(resultsfile,lotsizelabel,'Landscape Results','A43');
xlswrite(resultsfile,finedistances,'Landscape Results','B42');
xlswrite(resultsfile,meanlotsize,'Landscape Results','B43');
xlswrite(resultsfile,hpaclabel,'Landscape Results','A44');
xlswrite(resultsfile,meanhouseperacre,'Landscape Results','B44');
xlswrite(resultsfile,pllabel,'Landscape Results','A45');
xlswrite(resultsfile,meanplandint,'Landscape Results','B45');

xlswrite(resultsfile,distlabel_c,'Landscape Results','B47');
xlswrite(resultsfile,lotsizelabel,'Landscape Results','A49');
xlswrite(resultsfile,finedistances,'Landscape Results','B48');
xlswrite(resultsfile,meanlotsize_c,'Landscape Results','B49');
xlswrite(resultsfile,hpaclabel,'Landscape Results','A50');
xlswrite(resultsfile,meanhouseperacre_c,'Landscape Results','B50');
xlswrite(resultsfile,pllabel,'Landscape Results','A51');
xlswrite(resultsfile,meanplandint_c,'Landscape Results','B51');

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
xlswrite(resultsfile,lotlabel,'North Results','A12');
xlswrite(resultsfile,lottypeset,'North Results','B3');
xlswrite(resultsfile,lottypeset,'North Results','B12');
xlswrite(resultsfile,nrentstats(:,:,1),'North Results','C3');
xlswrite(resultsfile,nrentstats(:,:,2),'North Results','I3');
xlswrite(resultsfile,nrentstats(:,:,3),'North Results','O3');
xlswrite(resultsfile,nrentstats(:,:,4),'North Results','U3');
xlswrite(resultsfile,nrentstats(:,:,5),'North Results','AA3');
xlswrite(resultsfile,unzonedlotstats(:,:,1),'North Results','C12');
xlswrite(resultsfile,unzonedlotstats(:,:,2),'North Results','I12');
xlswrite(resultsfile,unzonedlotstats(:,:,3),'North Results','O12');
xlswrite(resultsfile,unzonedlotstats(:,:,4),'North Results','U12');
xlswrite(resultsfile,unzonedlotstats(:,:,5),'North Results','AA12');

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
xlswrite(resultsfile,lottypeset,'South Results','B12');
xlswrite(resultsfile,srentstats(:,:,1),'South Results','C3');
xlswrite(resultsfile,srentstats(:,:,2),'South Results','I3');
xlswrite(resultsfile,srentstats(:,:,3),'South Results','O3');
xlswrite(resultsfile,srentstats(:,:,4),'South Results','U3');
xlswrite(resultsfile,srentstats(:,:,5),'South Results','AA3');
xlswrite(resultsfile,zonedlotstats(:,:,1),'South Results','C12');
xlswrite(resultsfile,zonedlotstats(:,:,2),'South Results','I12');
xlswrite(resultsfile,zonedlotstats(:,:,3),'South Results','O12');
xlswrite(resultsfile,zonedlotstats(:,:,4),'South Results','U12');
xlswrite(resultsfile,zonedlotstats(:,:,5),'South Results','AA12');

xlswrite(resultsfile,timelabel,'Coastal Results','B1');
xlswrite(resultsfile,testtime(1),'Coastal Results','C1');
xlswrite(resultsfile,testtime(2),'Coastal Results','I1');
xlswrite(resultsfile,testtime(3),'Coastal Results','O1');
xlswrite(resultsfile,testtime(4),'Coastal Results','U1');
xlswrite(resultsfile,testtime(5),'Coastal Results','AA1');
xlswrite(resultsfile,lottypelabel,'Coastal Results','B2');
xlswrite(resultsfile,statlabel,'Coastal Results','C2');
xlswrite(resultsfile,statlabel,'Coastal Results','I2');
xlswrite(resultsfile,statlabel,'Coastal Results','O2');
xlswrite(resultsfile,statlabel,'Coastal Results','U2');
xlswrite(resultsfile,statlabel,'Coastal Results','AA2');
xlswrite(resultsfile,rentlabel,'Coastal Results','A3');
xlswrite(resultsfile,lotlabel,'Coastal Results','A12');
xlswrite(resultsfile,lottypeset,'Coastal Results','B3');
xlswrite(resultsfile,lottypeset,'Coastal Results','B12');
xlswrite(resultsfile,crentstats(:,:,1),'Coastal Results','C3');
xlswrite(resultsfile,crentstats(:,:,2),'Coastal Results','I3');
xlswrite(resultsfile,crentstats(:,:,3),'Coastal Results','O3');
xlswrite(resultsfile,crentstats(:,:,4),'Coastal Results','U3');
xlswrite(resultsfile,crentstats(:,:,5),'Coastal Results','AA3');
xlswrite(resultsfile,coastlotstats(:,:,1),'Coastal Results','C12');
xlswrite(resultsfile,coastlotstats(:,:,2),'Coastal Results','I12');
xlswrite(resultsfile,coastlotstats(:,:,3),'Coastal Results','O12');
xlswrite(resultsfile,coastlotstats(:,:,4),'Coastal Results','U12');
xlswrite(resultsfile,coastlotstats(:,:,5),'Coastal Results','AA12');

xlswrite(resultsfile,timelabel,'Midland Results','B1');
xlswrite(resultsfile,testtime(1),'Midland Results','C1');
xlswrite(resultsfile,testtime(2),'Midland Results','I1');
xlswrite(resultsfile,testtime(3),'Midland Results','O1');
xlswrite(resultsfile,testtime(4),'Midland Results','U1');
xlswrite(resultsfile,testtime(5),'Midland Results','AA1');
xlswrite(resultsfile,lottypelabel,'Midland Results','B2');
xlswrite(resultsfile,statlabel,'Midland Results','C2');
xlswrite(resultsfile,statlabel,'Midland Results','I2');
xlswrite(resultsfile,statlabel,'Midland Results','O2');
xlswrite(resultsfile,statlabel,'Midland Results','U2');
xlswrite(resultsfile,statlabel,'Midland Results','AA2');
xlswrite(resultsfile,rentlabel,'Midland Results','A3');
xlswrite(resultsfile,lotlabel,'Midland Results','A12');
xlswrite(resultsfile,lottypeset,'Midland Results','B3');
xlswrite(resultsfile,lottypeset,'Midland Results','B12');
xlswrite(resultsfile,midrentstats(:,:,1),'Midland Results','C3');
xlswrite(resultsfile,midrentstats(:,:,2),'Midland Results','I3');
xlswrite(resultsfile,midrentstats(:,:,3),'Midland Results','O3');
xlswrite(resultsfile,midrentstats(:,:,4),'Midland Results','U3');
xlswrite(resultsfile,midrentstats(:,:,5),'Midland Results','AA3');
xlswrite(resultsfile,middlelotstats(:,:,1),'Midland Results','C12');
xlswrite(resultsfile,middlelotstats(:,:,2),'Midland Results','I12');
xlswrite(resultsfile,middlelotstats(:,:,3),'Midland Results','O12');
xlswrite(resultsfile,middlelotstats(:,:,4),'Midland Results','U12');
xlswrite(resultsfile,middlelotstats(:,:,5),'Midland Results','AA12');

xlswrite(resultsfile,timelabel,'Inland Results','B1');
xlswrite(resultsfile,testtime(1),'Inland Results','C1');
xlswrite(resultsfile,testtime(2),'Inland Results','I1');
xlswrite(resultsfile,testtime(3),'Inland Results','O1');
xlswrite(resultsfile,testtime(4),'Inland Results','U1');
xlswrite(resultsfile,testtime(5),'Inland Results','AA1');
xlswrite(resultsfile,lottypelabel,'Inland Results','B2');
xlswrite(resultsfile,statlabel,'Inland Results','C2');
xlswrite(resultsfile,statlabel,'Inland Results','I2');
xlswrite(resultsfile,statlabel,'Inland Results','O2');
xlswrite(resultsfile,statlabel,'Inland Results','U2');
xlswrite(resultsfile,statlabel,'Inland Results','AA2');
xlswrite(resultsfile,rentlabel,'Inland Results','A3');
xlswrite(resultsfile,lotlabel,'Inland Results','A12');
xlswrite(resultsfile,lottypeset,'Inland Results','B3');
xlswrite(resultsfile,lottypeset,'Inland Results','B12');
xlswrite(resultsfile,inrentstats(:,:,1),'Inland Results','C3');
xlswrite(resultsfile,inrentstats(:,:,2),'Inland Results','I3');
xlswrite(resultsfile,inrentstats(:,:,3),'Inland Results','O3');
xlswrite(resultsfile,inrentstats(:,:,4),'Inland Results','U3');
xlswrite(resultsfile,inrentstats(:,:,5),'Inland Results','AA3');
xlswrite(resultsfile,inlandlotstats(:,:,1),'Inland Results','C12');
xlswrite(resultsfile,inlandlotstats(:,:,2),'Inland Results','I12');
xlswrite(resultsfile,inlandlotstats(:,:,3),'Inland Results','O12');
xlswrite(resultsfile,inlandlotstats(:,:,4),'Inland Results','U12');
xlswrite(resultsfile,inlandlotstats(:,:,5),'Inland Results','AA12');
% 
% obslabel={'Obs.'};
% avgmaplabel={'Avg Map'};
% pctdevlabel={'%dev'};
% threshlabel={'thresh'};
% xlswrite(resultsfile,obslabel,'AvgMap15','B25');
% xlswrite(resultsfile,avgmaplabel,'AvgMap15','C25');
% xlswrite(resultsfile,pctdevlabel,'AvgMap15','A26');
% xlswrite(resultsfile,threshlabel,'AvgMap15','A27');
% xlswrite(resultsfile,devstats(2,1),'AvgMap15','B26');
% xlswrite(resultsfile,testpctdev(2),'AvgMap15','C26');
% xlswrite(resultsfile,avgthresh(2),'AvgMap15','B27');
% 
% xlswrite(resultsfile,obslabel,'AvgMap20','B25');
% xlswrite(resultsfile,avgmaplabel,'AvgMap20','C25');
% xlswrite(resultsfile,pctdevlabel,'AvgMap20','A26');
% xlswrite(resultsfile,threshlabel,'AvgMap20','A27');
% xlswrite(resultsfile,devstats(3,1),'AvgMap20','B26');
% xlswrite(resultsfile,testpctdev(3),'AvgMap20','C26');
% xlswrite(resultsfile,avgthresh(3),'AvgMap20','B27');
% 
% xlswrite(resultsfile,obslabel,'AvgMap25','B25');
% xlswrite(resultsfile,avgmaplabel,'AvgMap25','C25');
% xlswrite(resultsfile,pctdevlabel,'AvgMap25','A26');
% xlswrite(resultsfile,threshlabel,'AvgMap25','A27');
% xlswrite(resultsfile,devstats(4,1),'AvgMap25','B26');
% xlswrite(resultsfile,testpctdev(4),'AvgMap25','C26');
% xlswrite(resultsfile,avgthresh(4),'AvgMap25','B27');
% 
% xlswrite(resultsfile,obslabel,'AvgMap30','B25');
% xlswrite(resultsfile,avgmaplabel,'AvgMap30','C25');
% xlswrite(resultsfile,pctdevlabel,'AvgMap30','A26');
% xlswrite(resultsfile,threshlabel,'AvgMap30','A27');
% xlswrite(resultsfile,devstats(5,1),'AvgMap30','B26');
% xlswrite(resultsfile,testpctdev(5),'AvgMap30','C26');
% xlswrite(resultsfile,avgthresh(5),'AvgMap30','B27');

xlswrite(resultsfile,timelabel,'Return Stats','B1');
xlswrite(resultsfile,(TSTART+1:TMAX),'Return Stats','C1');
xlswrite(resultsfile,lottypelabel,'Return Stats','B2');
xlswrite(resultsfile,retlabel,'Return Stats','A3');
xlswrite(resultsfile,sigmalabel,'Return Stats','A12');
xlswrite(resultsfile,lottypeset,'Return Stats','B3');
xlswrite(resultsfile,lottypeset,'Return Stats','B12');
xlswrite(resultsfile,meanlabel,'Return Stats','C2');
xlswrite(resultsfile,reshape(retstats(:,1,:),[HT 20]),'Return Stats','C3');
xlswrite(resultsfile,reshape(retstats(:,2,:),[HT 20]),'Return Stats','C12');

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
xlswrite(resultsfile,sigmalabel,'HtIncome Stats','A12');
xlswrite(resultsfile,lottypeset,'HtIncome Stats','B3');
xlswrite(resultsfile,lottypeset,'HtIncome Stats','B12');
xlswrite(resultsfile,meanlabel,'HtIncome Stats','C2');
xlswrite(resultsfile,reshape(htincomestats(:,1,:),[HT 20]),'HtIncome Stats','C3');
xlswrite(resultsfile,reshape(htincomestats(:,2,:),[HT 20]),'HtIncome Stats','C12');

xlswrite(resultsfile,timelabel,'BidShare Stats','B1');
xlswrite(resultsfile,(TSTART+1:TMAX),'BidShare Stats','C1');
xlswrite(resultsfile,lottypelabel,'BidShare Stats','B2');
xlswrite(resultsfile,bidlabel,'BidShare Stats','A3');
xlswrite(resultsfile,sigmalabel,'BidShare Stats','A12');
xlswrite(resultsfile,lottypeset,'BidShare Stats','B3');
xlswrite(resultsfile,lottypeset,'BidShare Stats','B12');
xlswrite(resultsfile,meanlabel,'BidShare Stats','C2');
xlswrite(resultsfile,reshape(bidsharestats(:,1,:),[HT 20]),'BidShare Stats','C3');
xlswrite(resultsfile,reshape(bidsharestats(:,2,:),[HT 20]),'BidShare Stats','C12');

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

xlswrite(resultsfile,ideallabel,'HouseOfferings','A14');
xlswrite(resultsfile,timelabel,'HouseOfferings','D15');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','A16');
xlswrite(resultsfile,lottypeset,'HouseOfferings','A17');
xlswrite(resultsfile,testtime,'HouseOfferings','B16');
xlswrite(resultsfile,avgidealset,'HouseOfferings','B17');

xlswrite(resultsfile,proflabel,'HouseOfferings','A26');
xlswrite(resultsfile,timelabel,'HouseOfferings','D27');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','A28');
xlswrite(resultsfile,lottypeset,'HouseOfferings','A29');
xlswrite(resultsfile,testtime,'HouseOfferings','B28');
xlswrite(resultsfile,avgprofset,'HouseOfferings','B29');

xlswrite(resultsfile,difflabel,'HouseOfferings','H14');
xlswrite(resultsfile,timelabel,'HouseOfferings','K15');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','H16');
xlswrite(resultsfile,lottypeset,'HouseOfferings','H17');
xlswrite(resultsfile,testtime,'HouseOfferings','I16');
xlswrite(resultsfile,diffset,'HouseOfferings','I17');

xlswrite(resultsfile,difflabel,'HouseOfferings','H26');
xlswrite(resultsfile,timelabel,'HouseOfferings','K27');
xlswrite(resultsfile,lottypelabel,'HouseOfferings','H28');
xlswrite(resultsfile,lottypeset,'HouseOfferings','H29');
xlswrite(resultsfile,testtime,'HouseOfferings','I28');
xlswrite(resultsfile,diffprofset,'HouseOfferings','I29');

