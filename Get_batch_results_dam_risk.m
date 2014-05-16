

% DAMAGE, lotchoice, farmsoldinfo, pctdev
MRUNS=25; 
NLENGTH=80;
NWIDTH=80;
load dist2coast
z=[0.5 1500; 0.5 2000; 0.5 2500; 2 1500; 2 2000; 2 2500; 5 1500; 5 2000; 5 2500];

coastvalue=linspace(1,100,MRUNS);
riskprcp=repmat(linspace(0,1,sqrt(MRUNS)),1,5);
dampct=reshape(repmat(linspace(0,1,sqrt(MRUNS)),5,1),1,MRUNS);
strmthrsh=linspace(0,100,MRUNS);

loc=zeros(NLENGTH,NWIDTH);
landprox=zeros(NLENGTH,NWIDTH);

landsaleinfo=cell(MRUNS,75); %[sale_price]
% totdam=cell(MRUNS,1);   %1x75
% housesaleinfo=cell(MRUNS,1);
% popinfo=cell(MRUNS,1); %1x75
% landsaleinfo=zeros(MRUNS,1); %[sale_price]
totdam=zeros(MRUNS,1);   %1x75
housesaleinfo=zeros(MRUNS,1);
popinfo=zeros(MRUNS,1); %1x75
avglandsales=zeros(MRUNS,1);
prox=1:75;
landsales=zeros(MRUNS,length(prox));
housepricetrans=zeros(1,NWIDTH);
fnames=dir;
fnamescell=struct2cell(fnames);
 %%% CHANGE FILE NAME %%%%
h=strncmp('results_CHALMS_Coast_batch_dam_risk',fnamescell(1,:),35);
% char(fnamescell{1,h})
fileorder=[1 12 22:28 2:11 13:21];
fnames=fnamescell(1,h);
for mr=1:MRUNS
    filename=fnames{fileorder(mr)};
    load(eval('filename'))
    
%     incometest=zeros(80,80);
%     riskmap=zeros(80,80);
%     inotvac=find(cat(1,lotchoice{:,4})==1);
%     for i=1:length(cat(1,Lottype{inotvac,1}))
%         incometest(Lottype{inotvac(i),2})=CONINFO{lotchoice{inotvac(i),5},1}*...
%             ones(length(Lottype{inotvac(i),2}),1);
%         riskmap(Lottype{inotvac(i),2})=CONINFO{lotchoice{inotvac(i),5},7}*...
%             ones(length(Lottype{inotvac(i),2}),1);
%     end
    
    farmerid=cat(2,LANDINFO{1,10});
    daminfo=cat(1,DAMAGE{:});
    daminfo=reshape(daminfo,NLENGTH,NWIDTH);
%     totdam(mr)=mat2cell(sum(daminfo(:,6:NWIDTH),1),1,75); %proximity=1:75
    totdam(mr)=sum(sum(daminfo(:,6:NWIDTH),1));
    
    numpop=reshape(LOTTYPE(:,30),NLENGTH,NWIDTH);
    iislot=(numpop~=0);
    ltcount=1./z(numpop(iislot),1);
    numpop(iislot)=ltcount;
%     popinfo(mr)=mat2cell(sum(numpop(:,6:NWIDTH),1),1,75);
    popinfo(mr)=sum(sum(numpop(:,6:NWIDTH),1));
    
    loc(lotlocate(:,2))=cat(1,lotchoice{lotlocate(:,1),7});
    for nl=6:NWIDTH
        inzero=(loc(:,nl)~=0);
        housepricetrans(nl)=median(loc(inzero,nl));
    end
%     housesaleinfo(mr)=mat2cell(housepricetrans(6:NWIDTH),1,75);
    housesaleinfo(mr)=median(housepricetrans(6:NWIDTH));
    
    for nf=1:length(farmsoldinfo(:,1))
        nfland=find(farmerid==farmsoldinfo(nf,1));
        uniqueid=unique(coastdist(nfland));
        for il=1:length(uniqueid)
            lsinfo=[cat(1,landsaleinfo{mr,uniqueid(il)}); ...
                farmsoldinfo(nf,4)*ones(length(find(coastdist(nfland)==uniqueid(il))),1)];
            landsaleinfo(mr,uniqueid(il))=mat2cell(lsinfo,length(lsinfo),1);
%             landsaleinfo(mr,uniqueid(il))=mat2cell(farmsoldinfo(nf,4)*...
%                 ones(length(find(coastdist(nfland)==uniqueid(il))),1), ...
%                 length(find(coastdist(nfland)==uniqueid(il))),1);
        end
    end
    for ip=1:length(prox)
        landsales(mr,ip)=mean(landsaleinfo{mr,ip});
    end
    ireal=~isnan(landsales(mr,:));
    avglandsales(mr)=mean(landsales(mr,ireal));
    
end
riskprcp=linspace(0,1,sqrt(MRUNS));
dampct=linspace(0,1,sqrt(MRUNS));
% Number of lots vs. prox
figure(1)
image(riskprcp,dampct,reshape(popinfo,5,5),'CDataMapping','scaled')
colorbar
axis xy
set(gcf,'Color','white')
xlabel('Risk Perception')
ylabel('Storm Severity')
title('Number of Residential Lots')

% Total Damages
figure(2)
image(riskprcp,dampct,reshape(totdam,5,5),'CDataMapping','scaled')
colorbar
axis xy
set(gcf,'Color','white')
xlabel('Risk Perception')
ylabel('Storm Severity')
title('Total Damages')

% Housing Prices
figure(3)
image(riskprcp,dampct,reshape(housesaleinfo,5,5),'CDataMapping','scaled')
colorbar
axis xy
set(gcf,'Color','white')
xlabel('Risk Perception')
ylabel('Storm Severity')
title('Median (Annualized) Housing Price')

% Land Prices
figure(4)
image(riskprcp,dampct,reshape(avglandsales,5,5),'CDataMapping','scaled')
colorbar
axis xy
set(gcf,'Color','white')
xlabel('Risk Perception')
ylabel('Storm Severity')
title('Mean Land Price')
