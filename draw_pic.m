function [pic3]=draw_pic(result)

pic3=[];
for sn=1:size(result,2)
    pic=repmat(result(:,sn),1,20);
    pic2=[];
    for i=1:size(result,1)
        pic2=[repmat(pic(i,:),3,1);pic2];
    end
    pic3=[pic3,pic2,zeros(size(result,1)*3,10)];
end
end
