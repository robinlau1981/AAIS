function [d,index]=mini_distance(X,Y)
[r,c]=size(Y);
for i=1:r
    if ~(X==Y(i,:))
       dis(i)=(X-Y(i,:))*(X-Y(i,:))';
    else
        dis(i)=inf;
    end
end
d=min(dis);
index=find(dis==d);