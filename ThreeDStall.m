function [Alpha2,CL3D,CD3D,AlphaZero]=ThreeDStall(RotorSPM,RotorR,WindV,rRatio,Chord_AirfoilPrep,InterpAlpha,InterpCl,InterpCd)

% ' Program to apply Du-Selig and Eggars 3-D stall corrections to 2-D airfoil data
% ' C Hansen, Windward Engineering, Dec 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Selig variables
RotorSPM;
RotorR;
WindV;
rRatio;
Chord_AirfoilPrep;
InterpAlpha;
InterpCl;
InterpCd;

A_constant=1;
B_constant=1;
D_constant=1;

CN2D=[];
CT2D=[];
CL3D=[]; 
CD3D=[];
CLP=[];
DelCL1=[];
DelCL2=[];
DelCD=[];
AlphaRad=[];

DToR=pi/180;

%'Get range of angles to be used for CL slope calc.
AlphaMinTrend=-4;
AlphaMaxTrend=4;

%Function UpdateCLSlope and find CL_slope 
X1=find(InterpAlpha == AlphaMinTrend);
X2=find(InterpAlpha == AlphaMaxTrend);
if abs(X1-X2)>40
    warning('You have requested more than the limit of 40 points for the CL Slope calc.')
end
Y_CL=InterpCl(X1:X2);
X_Alpha=InterpAlpha(X1:X2);
Slope= polyfit(X_Alpha, Y_CL, 1);
CLSlope= Slope(1);             %CL slope
Yintercept=Slope(2);            %intercept (y)
AlphaZero=-Yintercept/CLSlope;      %Alpha0
AlphaEnd=InterpAlpha(end,1);

%calculate and display some constants
cOverr=Chord_AirfoilPrep/RotorR/rRatio;
Lambda=RotorSPM*pi/30*RotorR/(WindV^2+(RotorSPM*pi/30*RotorR)^2)^0.5;
Expon = D_constant/Lambda/rRatio;
FL = 1/(CLSlope/DToR)*((1.6*cOverr/0.1267)*(A_constant-cOverr^Expon)/(B_constant+cOverr^Expon)-1);

%calculate 3D values
for i=1:length(InterpAlpha)
   
       AlphaRad(i)=DToR*InterpAlpha(i);      %' Alpha in radians
       CN2D(i)=InterpCl(i)*cos(AlphaRad(i))+InterpCd(i)*sin(AlphaRad(i));
       CT2D(i)=InterpCl(i)*sin(AlphaRad(i))-InterpCd(i)*cos(AlphaRad(i));
       
       CLP(i)=CLSlope*(InterpAlpha(i)-AlphaZero);
       
       DelCL1(i)=FL*(CLP(i)-InterpCl(i));
       
       if InterpAlpha(i)>AlphaEnd
          adj=((90-InterpAlpha(i))/(90-AlphaEnd))^2;
       else
          adj = 1;
       end
          
       CL3D(i)=InterpCl(i)+DelCL1(i)*adj;
       DelCL2(i)=CL3D(i)-InterpCl(i);
       
       DelCD(i)=DelCL2(i)*(sin(AlphaRad(i))-0.12*cos(AlphaRad(i)))/(cos(AlphaRad(i))+0.12*sin(AlphaRad(i)));
       CD3D(i)=InterpCd(i)+DelCD(i);
       
end
Alpha2=InterpAlpha;
CL3D=CL3D'; 
CD3D=CD3D';

