function [rv1,lv1]=auto_roi_mbr(cinemri)

[RV,LV]=FindLVRV(cinemri);

rv1=zeros(size(cinemri(:,:,15)));
lv1=zeros(size(cinemri(:,:,15)));

for i=-3:1:3
    for j=-3:1:3
        rv1(RV(1)+i+1,RV(2)+j+1)=1;
    end
end

for i=-3:1:3
    for j=-3:1:3
        lv1(LV(1)+i+1,LV(2)+j+1)=1;
    end
end


   
