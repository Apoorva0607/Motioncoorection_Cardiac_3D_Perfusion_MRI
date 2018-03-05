function murase_fitg(curves,slice)

% warning off all
% if nargin==0,
%     error('sliceNum argument needed for mpi_fitCurves');
% end

bldcurve=curves(1,:)';

nRegs=size(curves,1)-1;
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);

max_delay=round(size(curves,2))-1;

vare=zeros(size(curves,2),size(curves,1)-1);

num_par=3;
time=(0:(length(bldcurve)-1))/size(curves,2);


parfor ii=1:nRegs
        
    err=zeros(max_delay+1,1);    
    tissue=tisscurve(:,ii);
    
    for j=-max_delay:1:0
        
        t_delay=j;
        tmpbldcurve=bldcurve;
        tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
        tmpbldcurve(1:1-fix(t_delay))=0;
        try
            [~,fit,~] = muraseFit(time,tmpbldcurve,tissue,num_par);        
        catch
            fit=bldcurve;
        end
        err(j+max_delay+1)=sum((tissue-fit).^2);        
    end
    
    pos=find(err==min(err));
    t_delay=(pos)-(max_delay+1);
    
    tmpbldcurve=bldcurve;
    tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
    tmpbldcurve(1:1-fix(t_delay))=0;

    tissue=tisscurve(:,ii);    
    [~,fit,~] = muraseFit(time,tmpbldcurve,tissue,num_par);
    
    vare(:,ii)=fit;

end


    temp=strcat('Output1/',int2str(slice),'/vare.mat');
save(temp,'vare');
return;

