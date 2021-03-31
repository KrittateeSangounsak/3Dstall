function [StallAngle,NegStall,CNSlope,CN1,CN2,CDMin]=Dynstall(Alpha2,CL2,CD2,nTable2)
  
%Dynamic stall parameters
CNSlope=[];
AlphaCDMin=[];
CNMax=[];

%data from TableExtrap
CN=[];
CT=[];
nTable1=[];
Count=0;
AlphaMinTrend=-2;
AlphaMaxTrend=8;

DToR=pi/180;

%calculate CN, CT
for i=1:nTable2
    CN(i)=CL2(i)*cos(DToR*Alpha2(i))+CD2(i)*sin(DToR*Alpha2(i));
    CT(i)=CL2(i)*sin(DToR*Alpha2(i))-CD2(i)*cos(DToR*Alpha2(i));
        
end
CN=CN';
CT=CT';

%Function UpdateCNSlope and find CN_slope 
X1=find(Alpha2 == AlphaMinTrend);
X2=find(Alpha2 == AlphaMaxTrend);
if abs(X1-X2)>40
    warning('You have requested more than the limit of 40 points for the CN Slope calc.')
end
Y_CN=CN(X1:X2);
X_Alpha=Alpha2(X1:X2);
Slope= polyfit(X_Alpha, Y_CN, 1);
CNSlope= Slope(1)*(180/pi);

%get stall angle 
[CNMax,Index]=max(CN);
CNMax;
StallAngle=Alpha2(Index);
NegStall=-Alpha2(Index);

%get Zero Cn angle of attack (deg)
ZeroCn=-(Slope(2)/CNSlope)*(180/pi);

%Cn extrapolated to value at positive stall angle of attack
CN1=(CNSlope*pi/180*(StallAngle-ZeroCn));

%Cn at stall value for negative angle of attack
CN2 = interp1(Alpha2,CN,NegStall);

%get angle of attack for CDmin 
%(search only in range -20 to +20?)
Z1=find(Alpha2 == -20);
Z2=find(Alpha2 == 20);
CDmin=CD2(Z1:Z2);
CDMin = 100;
for i=Z1:Z2
    if CD2(i) <= CDMin
    CDMin = CD2(i);
    AlphaCDMin= Alpha2(i);
    end 
end
CDMin;
AlphaCDMin;








