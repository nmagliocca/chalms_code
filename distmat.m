NLENGTH=80;
NWIDTH=80;
NCELLS=NLENGTH*NWIDTH;
ddist2hznnei=zeros(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
ddist2vrtnei=zeros(NLENGTH,NWIDTH);
distmat=cell(NCELLS,1);
for tl=1:NCELLS
    %%% Distance calc needs to be revised when applied to irregular grid
    inddist2dev=10000*ones(NLENGTH,NWIDTH);
    [vacrow,vaccol]=ind2sub([NLENGTH NWIDTH],tl);
    
    for col=1:NWIDTH
        ddist2hznnei(1:NLENGTH,col)=abs(col-vaccol).*...
            ones(NLENGTH,1);
        for row=1:NWIDTH
            ddist2vrtnei(row,1:NWIDTH)=abs(row-vacrow).*...
                ones(1,NLENGTH);
            inddist2dev(row,col)=min(sqrt(ddist2hznnei(row,col)^2+...
                ddist2vrtnei(row,col)^2),inddist2dev(row,col));
        end
    end
    distmat(tl)=mat2cell(reshape(inddist2dev,NCELLS,1),NCELLS,1);
end
