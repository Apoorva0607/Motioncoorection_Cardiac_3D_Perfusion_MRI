function cine_curv_fit(vart,noi,Y,X,slice)

temp=strcat('Output1/',int2str(slice),'/init_bld.mat');
load(temp);
temp=strcat('Output1/',int2str(slice),'/init_tiss.mat');
load(temp);
temp=strcat('Output1/',int2str(slice),'/tiss_avgSpatial.mat');
load(temp);

if(vart==1)
temp=strcat('Output1/',int2str(slice),'/varg'); % fmincon
load(temp);
end

if(vart==0)
temp=strcat('Output1/',int2str(slice),'/vare');; % murase
load(temp)
end

if(vart==1)
est_cur=varg; 
else
    est_cur=vare;
end

nRegs=size(est_cur,2);

orig_tisscurve=zeros(size(est_cur));

for i=1:nRegs
 orig_tisscurve(:,i)=(est_cur(:,i)*init_tiss(i))/(tiss_avgSpatial)+init_tiss(i);
end

temp=strcat('Output1/',int2str(slice),'/cinemri1.mat');
load(temp);

cinemri=cinemri1;

img=cinemri;
new_tisscurves=orig_tisscurve';

cinemri_curv_fit=img;
nX=length(X);
not_fit=200*ones(1,size(new_tisscurves,2));
counter=0;

for j = 1 : nX
    a=find(new_tisscurves(j,:)<10000);
    if(a)
      cinemri_curv_fit(Y(j), X(j), :)= new_tisscurves(j,:);
  else
      counter=counter+1;
      cinemri_curv_fit(Y(j), X(j), :)=not_fit;
    end
end
     
temp=strcat('Output1/',int2str(slice),'/cinemri_curv_fit',int2str(noi),'.mat');

save(temp,'cinemri_curv_fit');

