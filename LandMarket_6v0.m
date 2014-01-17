%%%%%%%%%%%%%% ABM Land Market %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdepart=zeros([],1);

for tt=t
    if isempty(iNfarmers)==1 || landdemand(tt) == 0
        ifarmtrans=[];
%         Plandproj(iNfarmers,tt)=mean([wtaland(iNfarmers,tt) wtpland(iNfarmers,tt)],2);
        continue
    else
        ipool=find(wtaland(iNfarmers,tt)<=wtpland(iNfarmers,tt));
        ioutpool=find(wtaland(iNfarmers,tt)>=wtpland(iNfarmers,tt));
    end

    if isempty(ipool)==1
        ifarmtrans=[];
%         Plandproj(iNfarmers,tt)=mean([wtaland(iNfarmers,tt) wtpland(iNfarmers,tt)],2);
        Plandproj(iNfarmers,tt)=0.75*wtaland(iNfarmers,tt)+0.25*wtpland(iNfarmers,tt);
        Paskland(iNfarmers,tt)=wtaland(iNfarmers,tt);
        Pdevbid(iNfarmers,tt)=wtpland(iNfarmers,tt);
        continue
    end
%     epsilon(tt)=zeta*(nproj(ibest)-sum(Farmstats(iNfarmers(ipool),2,t)))/...
%         (nproj(ibest)+sum(Farmstats(iNfarmers(ipool),2,t)));
    epsilon(tt)=zeta*(landdemand(tt)-sum(Farmstats(iNfarmers(ipool),2,tt)))/...
        (landdemand(tt)+sum(Farmstats(iNfarmers(ipool),2,tt)));

    Paskland(iNfarmers(ipool),tt)=max(wtaland(iNfarmers(ipool),tt).*...
        (1+epsilon(tt)),wtaland(iNfarmers(ipool),tt));
    Paskland(iNfarmers(ioutpool),tt)=wtaland(iNfarmers(ioutpool),tt);
    Plandproj(iNfarmers(ioutpool),tt)=wtpland(iNfarmers(ioutpool),tt);
    Pdevbid(iNfarmers(ipool),tt)=min(wtpland(iNfarmers(ipool),tt)*...
        (1+epsilon(tt)),wtpland(iNfarmers(ipool),tt));
    
    
    %%% Farmer Selection Criteria %%%
    isell=find(Paskland(iNfarmers(ipool),tt)<= Pdevbid(iNfarmers(ipool),tt));
    if isempty(isell)==1
        ifarmtrans=[];
        Plandproj(iNfarmers(ipool),tt)=Pdevbid(iNfarmers(ipool),tt);
%         Plandproj(iNfarmers(ipool),tt)=0.75*Paskland(iNfarmers...
%             (ipool),tt)+0.25*Pdevbid(iNfarmers(ipool),tt);
        continue
    else
        transprice(iNfarmers(ipool(isell)),tt)=mean([Paskland(iNfarmers...
            (ipool(isell)),tt) Pdevbid(iNfarmers(ipool(isell)),tt)],2);
%         transprice(iNfarmers(ipool(isell)),tt)=Paskland(iNfarmers(ipool(isell)),tt);
        
        for nnf=1:length(isell)
            ifarmcalc=find(AGLAYER == iNfarmers(ipool(isell(nnf))));
            maxrent(iNfarmers(ipool(isell(nnf))),tt)=sum(maxRENTPROJ(ifarmcalc));
%             maxreturn(iNfarmers(ipool(isell(nnf))),tt)=maxrent(iNfarmers(ipool(isell(nnf))),tt)...
%                 -(length(ifarmcalc)*(transprice(iNfarmers(ipool(isell(nnf))),tt)*discount)+...
%                 sum(CCOST(ifarmcalc)./ZZ(ifarmcalc)));
            maxreturn(iNfarmers(ipool(isell(nnf))),tt)=sum(maxret(ifarmcalc))-...
                discount*transprice(iNfarmers(ipool(isell(nnf))),tt).*...
                Farmstats(iNfarmers(ipool(isell(nnf))),2,1);
        end
        sellinfo=[maxreturn(iNfarmers(ipool(isell)),tt) ipool(isell) ...
            Farmstats(iNfarmers(ipool(isell)),2,1)];
        lowprice=sortrows(sellinfo,-1);
%         lowprice=sortrows(sellinfo,1);
        landsupply=cumsum(lowprice(:,3));
        ilandbuy=find(landsupply >= landdemand(tt),1,'first');
%         ifarmer=lowprice(:,2);
%         keyboard
        if isempty(ilandbuy)==0
            ibuy=1:ilandbuy;
        elseif isempty(ilandbuy)==1
            ibuy=1:length(isell);
        end
    end
   
    fdepart(length(fdepart)+1:length(fdepart)+length(ibuy),1)=...
        iNfarmers(lowprice(ibuy,2));
    for fd=1:length(fdepart)
        PLAND(AGLAYER==fdepart(fd))=transprice(fdepart(fd),tt);
    end
    sellrecord(fdepart,1)=tt;
    buyrecord(fdepart,1)=transprice(iNfarmers(lowprice(ibuy,2)),tt);
    iselltime=(sellrecord>0);

    ifarmtrans=fdepart;
end