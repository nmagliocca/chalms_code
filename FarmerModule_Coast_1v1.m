%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Price Projection models
learnlandproj = zeros(Nfarmers,NUMMODEL);
learnlanderror = zeros(Nfarmers,NUMMODEL);
learnlandbestSAVE=zeros(Nfarmers,TMAX);
ilearnlandbestSAVE=zeros(Nfarmers,TMAX);
learnlandprojSAVE=zeros(Nfarmers,NUMMODEL);
learnlandmodelSAVE=zeros(Nfarmers,TMAX);
learnlandprojbestSAVE=zeros(Nfarmers,NUMMODEL);
learnlandmodelbestSAVE=zeros(Nfarmers,TMAX);
learnlandprojbestSAVE=zeros(Nfarmers,TMAX);
learnwtaland=zeros(Nfarmers,TMAX);
learnDELTA=0.95;

projdiff=zeros(Nfarmers,TMAX);
pctdiff=zeros(Nfarmers,TMAX);
pricevar=zeros(Nfarmers,TMAX);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%    FARMERS LEARNING MODULE    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iinitialcity=(BASELAYER==1 & SCAPE==1);
firstsold=unique(AGLAYER(iinitialcity));
ifarmtrans=firstsold;       %co-opted from land market
meanPLAND(1)=mean(PLAND(iinitialcity));
FARMPROJ=meanPLAND(1)-travelcost;
for rf=1:length(iNfarmers)
    ifarmcalc=find(AGLAYER == iNfarmers(rf));
    sumret=FARMPROJ(ifarmcalc);
    wtpland(iNfarmers(rf),1)=sum(sumret)/length(ifarmcalc);
end
learnwtaland(iNfarmers,1)=wtaland(iNfarmers,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
learnwtaland(iNfarmers,2)=learnwtaland(iNfarmers,1);
successflag=1;
successcount=0;
successmark=zeros(1,[]);
tlearn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    Farmer Price Prediction    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while successcount <= 2
    tlearn=tlearn+1;
    %%% Linear, Random price signal
    PLAND(iinitialcity)=PLAND(iinitialcity)+200*rand(1);

    meanPLAND(tlearn)=mean(PLAND(iurb));
    FARMPROJ=meanPLAND(tlearn)-travelcost;
    for rf=1:length(iNfarmers)
        ifarmcalc=find(AGLAYER == iNfarmers(rf));
        sumret=FARMPROJ(ifarmcalc);
        wtpland(iNfarmers(rf),tlearn)=sum(sumret)/length(ifarmcalc);
    end
    
    %%%%% LAND MARKET %%%%%%%%

    Paskland(iNfarmers,tlearn)=learnwtaland(iNfarmers,tlearn);
    Plandproj(iNfarmers,tlearn)=mean([learnwtaland(iNfarmers,tlearn) ...
        wtpland(iNfarmers,tlearn)],2);

    % Calculate Distances to Transactions
    ntrans=length(ifarmtrans);   %number of land market trans actions, farms sold. In this module, for farmer learning, initial city farms sold held constatnt
    indtransdist=zeros(NLENGTH,NWIDTH,ntrans);
    transdist=zeros(NLENGTH,NWIDTH);
    Planddistproj=zeros(length(iNfarmers),ntrans);

    if ntrans > 0
        for nt=1:ntrans
            tdist2hznnei=10000*ones(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
            tdist2vrtnei=10000*ones(NLENGTH,NWIDTH);
            transland=find(AGLAYER==ifarmtrans(nt));
            for tl=1:length(transland)
                [itlrow,itlcol]=ind2sub([NLENGTH,NWIDTH],transland(tl));
                for col=1:NWIDTH
                    tdist2hznnei(1:NLENGTH,col)=min(abs(col-itlcol).*...
                        ones(NLENGTH,1),tdist2hznnei(1:NLENGTH,col));
                end
                for row=1:NWIDTH
                    tdist2vrtnei(row,1:NWIDTH)=min(abs(row-itlrow).*...
                        ones(1,NLENGTH),tdist2vrtnei(row,1:NWIDTH));
                end
                for col=1:NWIDTH
                    for row=1:NLENGTH
                        indtransdist(row,col,nt)=sqrt(tdist2hznnei(row,col)^2+tdist2vrtnei(row,col)^2);
                    end
                end
            end
            iNfarmers=unique(AGLAYER(iscape));
            %Calculate gradient coefficients using genetic algorithm
            subtransdist=indtransdist(:,:,nt);
%             distcoeff(ifarmtrans,:)=0;
            for nf=1:length(iNfarmers)
                rc=zeros([],1);
                avgtransdist=max(mean(subtransdist(AGLAYER==iNfarmers(nf))),1);
                coeffmark=(Paskland(iNfarmers(nf),tlearn)-mean(PLAND(transland)))/...
                    avgtransdist;
                fitness(iNfarmers(nf),:,tlearn)=fitness(iNfarmers(nf),:,tlearn-1)+...
                    abs(distcoeff(iNfarmers(nf),:)-coeffmark);
                fitsort=sort(fitness(iNfarmers(nf),:,tlearn),'ascend');
                stratcount=1;
                for x=1:nextgen
                    numstrat=find(fitness(iNfarmers(nf),:,tlearn)==fitsort(x));
                    if length(numstrat) > 1
                        isurvive(x,1)=numstrat(stratcount);
                        stratcount=stratcount+1;
                    else
                        isurvive(x,1)=numstrat(stratcount);
                    end
                    if stratcount > length(numstrat)
                        stratcount=1;
                    end
                end
                distcoeff(iNfarmers(nf),1:nextgen)=distcoeff(iNfarmers(nf),isurvive);
                for xx=1:nextgen
                    irc=xx*recombo/nextgen;
                    rc(length(rc)+1:irc,1)=distcoeff(iNfarmers(nf),xx)+rand(1,4);
                end
                distcoeff(iNfarmers(nf),nextgen+1:nextgen+recombo)=rc';
                distcoeff(iNfarmers(nf),nextgen+recombo+1:NUMMODELDIST)=...
                    mincoeff+(maxcoeff-mincoeff)*rand(1,NUMMODELDIST-...
                    (nextgen+recombo));
            end
            Planddistproj(iNfarmers,nt)=distcoeff(iNfarmers,1).*avgtransdist+...
                mean(PLAND(transland));
        end
        for nnw=1:NWIDTH
            for nnl=1:NLENGTH
                transdist(nnl,nnw)=min(indtransdist(nnl,nnw,1:ntrans));
            end
        end
        Planddist(iNfarmers,tlearn)=mean(Planddistproj(iNfarmers,:),2);
    else
        Planddist(iNfarmers,tlearn)=mean([learnwtaland(iNfarmers,tlearn) wtpland(iNfarmers,tlearn)],2);
    end
    %%% Land Price projections for tlearn+1
    for nf=1:length(iNfarmers)

        ilandclass1=find(landmodel(iNfarmers(nf),:)==1);
        ilandclass2=find(landmodel(iNfarmers(nf),:)==2);
        ilandclass3=find(landmodel(iNfarmers(nf),:)==3);
        ilandclass4=find(landmodel(iNfarmers(nf),:)==4);
        ilandclass5=find(landmodel(iNfarmers(nf),:)==5);
        ilandclass6=find(landmodel(iNfarmers(nf),:)==6);
        for i = 1:FARMNUMCLASS
            if i == 1
                % mimic models
                learnlandproj(iNfarmers(nf),ilandclass1) = Plandproj(iNfarmers(nf),tlearn)+(1-aa...
                    (iNfarmers(nf),ilandclass1)).*(0.5*Plandproj(iNfarmers(nf),tlearn)-...
                    (Plandproj(iNfarmers(nf),tlearn)-Plandproj(iNfarmers(nf),tlearn-1)));
            elseif i == 2
                % mean model
                for jl = 1:length(ilandclass2)
                    learnlandproj(iNfarmers(nf),ilandclass2(jl)) = mean(Plandproj(iNfarmers(nf),...
                        tlearn:-1:max((tlearn-aa(iNfarmers(nf),ilandclass2(jl))),1)));
                end
                
            elseif i == 3
                %cycle model
                learnlandproj(iNfarmers(nf),ilandclass3) = Plandproj(iNfarmers(nf),...
                    max(tlearn-round(aa(iNfarmers(nf),ilandclass3)),1));
            elseif i == 4
                % projection model
                for jl = 1:length(ilandclass4)
                    %Nonlinear Forecast
                    warning off
                    indata=Plandproj(iNfarmers(nf),tlearn-min(tlearn-1,(1+aa(iNfarmers(nf),ilandclass4(jl)))):tlearn);
%                     pcoef=polyfit(1:length(indata),indata,2);
%                     pline=pcoef(1).*(1:length(indata)+1).^2+pcoef(2).*(1:...
%                         length(indata)+1)+pcoef(3);
                    pcoef=polyfit(1:length(indata),indata,1);
                    pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
                    learnlandproj(iNfarmers(nf),ilandclass4(jl))=pline(length(pline));
                    warning on
                end
                
            elseif i == 5
                % rescale model
                learnlandproj(iNfarmers(nf),ilandclass5) = aa(iNfarmers(nf),...
                    ilandclass5)*Plandproj(iNfarmers(nf),tlearn);
            elseif i == 6
                % local(1) or regional(0) trends
                ilandreg=(aa(iNfarmers(nf),ilandclass6)==0);
                ilandlocal=(aa(iNfarmers(nf),ilandclass6)==1);
                % Local: just spatially discounted
                learnlandproj(iNfarmers(nf),ilandclass6(ilandlocal)) = Planddist...
                    (iNfarmers(nf),tlearn);
                if isempty(iNfarmers)==1
                    break
                end
                % Regional: density-dependent
                learnlandproj(iNfarmers(nf),ilandclass6(ilandreg)) = Planddist...
                    (iNfarmers(nf),tlearn).*(1+1/length(iNfarmers));
            end
        end

        if tlearn > 2
            learnlanderror(iNfarmers(nf),:) = (1-learnDELTA)*learnlanderror(iNfarmers(nf),:)+...
                learnDELTA*abs(Plandproj(iNfarmers(nf),tlearn)-learnlandprojSAVE(iNfarmers(nf),:));
            % Use model that predicted this period's price best
            [landbest ilandbest] = min(learnlanderror(iNfarmers(nf),:),[],2);
        else
            [landbest ilandbest] = min(learnlanderror(iNfarmers(nf),:),[],2);
        end

        learnwtaland(iNfarmers(nf),tlearn+1)=max(learnlandproj(iNfarmers(nf),ilandbest),...
            Farmstats(iNfarmers(nf),5,1));

        learnlandprojSAVE(iNfarmers(nf),:)=learnlandproj(iNfarmers(nf),:);
        learnlandbestSAVE(iNfarmers(nf),tlearn) = landbest;
        ilearnlandbestSAVE(iNfarmers(nf),tlearn) = ilandbest;
        learnlandprojbestSAVE(iNfarmers(nf),tlearn+1) = learnlandproj(iNfarmers(nf),ilandbest);
        learnlandmodelbestSAVE(iNfarmers(nf),tlearn) = landmodel(iNfarmers(nf),ilandbest);
    end

    %%%% End tlearn loop

    %         figure(1)
    %         hist(landmodelbestSAVE(:,tlearn),1:6)

    projdiff(:,tlearn)=Plandproj(:,tlearn)-learnlandprojbestSAVE(:,tlearn);
    pctdiff(:,tlearn)=abs(projdiff(:,tlearn))./Plandproj(:,tlearn);
    pricevar(:,tlearn)=abs(diff(Plandproj(:,tlearn-1:tlearn),1,2)./Plandproj(:,tlearn));
    if tlearn > TSTART+1
        %             successflag=(pctdiff(:,tlearn) >= pricevar(:,tlearn));
        %             successflag=(sum(pctdiff(:,tlearn)) >= sum(pricevar(:,tlearn)))+...
        %                 (length(find(pctdiff(:,tlearn) <= pricevar(:,tlearn))==1) <= 0.9*Nfarmers);
        successflag=(length(find(pctdiff(:,tlearn) < 0.10)==1) <= 0.9*Nfarmers);
        if successflag==1
            successmark(length(successmark)+1)=0;
            successcount=0;
        elseif successflag==0
            successmark(length(successmark)+1)=1;
            successcount=successcount+1;
        end
    end
end

% Feed model success information to TSTART

fitnessSAVE=fitness(:,:,tlearn);
