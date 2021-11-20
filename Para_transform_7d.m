function Yt=Para_transform_7d(Y)   
Yt=Y;
    Yt(2)=exp(Y(2));
    Yt(3)=exp(Y(3));
    Yt(4)=sqrt(Y(4)^2+Y(5)^2);
    Yt(5)=mod(atan(Y(5)/Y(4))+pi*(Y(4)<0),2*pi);
    Yt(6)=mod(Y(6)-Yt(5),2*pi);
    Yt(7)=exp(Y(7));