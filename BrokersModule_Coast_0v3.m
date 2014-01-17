%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    Broker's Learning Module    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Price Projection models
learnbproj=zeros(HT,NUMMODEL);
learnbrokerproj = zeros(HT,NUMMODEL,Nbrokers);
learnbrokererror = zeros(HT,NUMMODEL,Nbrokers,TMAX);
learnabserror = zeros(HT,NUMMODEL,Nbrokers,TMAX);
laernbrokerbestSAVE=zeros(Nbrokers,HT,TMAX);
ilearnbrokerbestSAVE=zeros(Nbrokers,HT,TMAX);
learnbrokerprojSAVE=zeros(Nbrokers,HT,TMAX);
learnbrokermodelSAVE=zeros(Nbrokers,HT,TMAX);
learnDELTA=0.5;
difflearnerror=zeros(HT,NUMMODEL,Nbrokers,TMAX);
learnbestdiffSAVE=zeros(Nbrokers,HT,TMAX);
learnbestabsSAVE=zeros(Nbrokers,HT,TMAX); 
projdiff=zeros(HT,Nbrokers,TMAX);
pctdiff=zeros(HT,Nbrokers,TMAX);
pricevar=zeros(HT,Nbrokers,TMAX);

% Input sample rents
testrents=10000*[0.6721    0.6999    0.8438    0.9093    0.9526
    0.8860    0.9349    1.1483    1.2432    1.2900
    1.2645    1.3129    1.5328    1.6199    1.6786
    1.9035    1.9304    1.9204    1.9743    2.0447
    2.0947    2.1348    2.1311    2.1762    2.2347
    2.4736    2.5188    2.5519    2.6037    2.6382
    02.35     2.37      2.42      2.46      2.5
    2.55      2.6       2.63      2.67      2.7
    2.83      2.89      2.94      2.96      3.0];
testtime=[10 15 20 25 30];
pcoeffs=zeros(HT,2);
for lt=1:HT
    pcoeffs(lt,:)=polyfit(testtime,testrents(lt,:),1);
end

% rentinfo=pcoeffs(1:12,1)*(1:TMAX)+pcoeffs(1:12,2)*ones(1,TMAX)+(500*randn(12,TMAX));
rentinfo=pcoeffs(1:HT,1)*(1:TMAX)+(pcoeffs(1:HT,2)*...
    ones(1,TMAX).*(0.85+(1.15-0.85)*rand(HT,TMAX)));

successflag=1;
successcount=0;
successmark=zeros(1,[]);
tlearn=TSTART;

while successcount <= 2
    tlearn=tlearn+1;

    for nb=1:Nbrokers
        
        bb=bbfull(:,:,nb);
        for lt=1:HT
            for i = 1:BROKERNUMCLASS
                %         strb = sprintf('brokerclass%d = find(brokermodel(:,:,nb) == %d);',i,i);
                strb = sprintf('brokerclass%d = find(brokermodel(lt,:,nb) == %d);',i,i);
                eval(strb);
            end

            for i = 1:BROKERNUMCLASS
                if i == 1
                    % mimic models
                    learnbproj(lt,brokerclass1) = rentinfo(lt,tlearn)+(1-bb...
                        (lt,brokerclass1)).*(0.5*rentinfo(lt,tlearn)-...
                        (rentinfo(lt,tlearn)-rentinfo(lt,tlearn-1)));
                elseif i == 2
                    % mean model
                    for jl = 1:length(brokerclass2)
                        learnbproj(lt,brokerclass2(jl)) = mean(rentinfo(lt,...
                            tlearn:-1:(tlearn-bb(lt,brokerclass2(jl)))));
                    end
                elseif i == 3
                    %cycle model
                    learnbproj(lt,brokerclass3) = rentinfo(lt,tlearn-...
                        round(max(1,bb(lt,brokerclass3))));
                elseif i == 4
                    % projection model
                    for jl = 1:length(brokerclass4)
                        %Nonlinear Forecast
                        indata=rentinfo(lt,tlearn-(1+bb(lt,brokerclass4(jl))):tlearn);
                        subindata=reshape(indata,1,length(indata));
                        pcoef=polyfit(1:length(indata),subindata,1);
                        pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
%                         pcoef=polyfit(1:length(indata),subindata,2);
%                         pline=pcoef(1).*(1:length(indata)+1).^2+pcoef(2).*(1:...
%                             length(indata)+1)+pcoef(3);
                        learnbproj(lt,brokerclass4(jl))=pline(length(pline));
                    end
                elseif i == 5
                    % rescale model
                    learnbproj(lt,brokerclass5) = bb(lt,brokerclass5)*rentinfo(lt,tlearn);
                elseif i == 6
                    [brows bcols]=ind2sub([nbrokerlong nbrokerwide],nb);
                    brnei=(bcols+1)*nbrokerlong-(nbrokerwide-brows);
                    blnei=(bcols-1)*nbrokerlong-(nbrokerwide-brows);
                    bupnei=bcols*nbrokerlong-(nbrokerwide-(brows-1));
                    bdnnei=bcols*nbrokerlong-(nbrokerwide-(brows+1));
                    ibnei=[brnei blnei bupnei bdnnei];
                    realbnei=find(minibmap==brnei | minibmap==blnei | ...
                        minibmap==bupnei | minibmap==bdnnei);
                    learnbproj(lt,brokerclass6) = (bb(lt,brokerclass6)+...
                        rand(1,length(brokerclass6)))*mean(rentinfo(lt,tlearn));
                end
            end

            learnbrokererror(lt,:,nb,tlearn) = (1-learnDELTA)*learnbrokererror(lt,:,nb,tlearn-1)+...
                learnDELTA*abs(rentinfo(lt,tlearn)-learnbproj(lt,:));
            learnabserror(lt,:,nb,tlearn)=rentinfo(lt,tlearn)-learnbproj(lt,:);
            [brokerbest ibrokerbest] = min(learnbrokererror(lt,:,nb,tlearn),[],2);
            if tlearn > TSTART+1
                difflearnerror(lt,:,nb,tlearn)=learnbrokererror(lt,:,nb,tlearn)-...
                    learnbrokererror(lt,:,nb,tlearn-1);
                learnbestdiffSAVE(nb,lt,tlearn)=difflearnerror(lt,ibrokerbest,nb,tlearn);
                learnbestabsSAVE(nb,lt,tlearn)=learnabserror(lt,ibrokerbest,nb,tlearn);
            else
                difflearnerror(lt,:,nb,tlearn)=0;
            end
            learnbrokerbestSAVE(nb,lt,tlearn) = brokerbest';
            ilearnbrokerbestSAVE(nb,lt,tlearn) = ibrokerbest';
            learnbrokerprojSAVE(nb,lt,tlearn) = learnbproj(lt,ibrokerbest);
            learnbrokermodelSAVE(nb,lt,tlearn) = brokermodel(lt,ibrokerbest,nb);

        end
        brokerproj(:,:,nb)=learnbproj;

        projdiff(:,nb,tlearn)=rentinfo(:,tlearn)-learnbrokerprojSAVE(nb,:,tlearn)';
        pctdiff(:,nb,tlearn)=abs(projdiff(:,nb,tlearn))./rentinfo(:,tlearn);
        pricevar(:,nb,tlearn)=abs(diff(rentinfo(:,tlearn-1:tlearn),1,2)./...
            rentinfo(:,tlearn));
    end
    if tlearn > TSTART+TSTART+1
        meansuccess=mean(pctdiff(:,:,tlearn));
        successflag=(length(find(meansuccess < 0.10)==1) <= 0.9*Nbrokers);
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

