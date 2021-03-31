function [Alpha2,CL2,CD2,CM2,CM0,nTable2]=Extrapolation(Alpha1,CL1,CD1,CM1)

nTable1=length(Alpha1);

DToR=pi/180;

%Function Viterna 
AspectRatio=17;
CDMax=min(2.01,1.11+0.018*AspectRatio);
CLRef=CL1(end,1);
CDRef=CD1(end,1);          %CL at upper matching poin
CMCoef=[];                 %Coefficient used in CM calcs.
CM0=[];
VAlphaLo=Alpha1(1,1);      %lower matching point
VAlphaHi=Alpha1(end,1);    %upper matching point
Alpha2=[];                 %new table
CL2=[];
CD2=[];
CM2=[];

%Get Viterna coefficients
SAlpha=sin(VAlphaHi*DToR);
CAlpha=cos(VAlphaHi*DToR);
A2 = (CLRef-CDMax/2*sin(2*VAlphaHi*DToR))*SAlpha/CAlpha^2;
B2 = (CDRef-CDMax*SAlpha^2)/CAlpha;

%Function CMCoeff
CLLo=CL1(1,1);
CDLo=CD1(1,1);
CLHi=CL1(end,1);
CDHi=CD1(end,1);
CMHi=CM1(end,1);
FoundZeroLift=0;
% Get CM at angle of zero lift (CM0)
for i=1:nTable1-1
    if abs(Alpha1(i))<20 && CL1(i)<=0 && CL1(i+1)>= 0
       p = -CL1(i)/(CL1(i + 1)-CL1(i));     %interpolation parameter
       CM0 = CM1(i)+ p*(CM1(i + 1)-CM1(i));
       FoundZeroLift = 1;
       break
    end
end

%zero lift not in range, use first two points
if FoundZeroLift==0 
        p = -CL1(1)/(CL1(2)-CL1(1));        %interpolation parameter
      CM0 = CM1(1)+ p*(CM1(2)-CM1(1));
end

XM =(-CMHi+CM0)/(CLHi*cos(VAlphaHi*DToR)+CDHi*sin(VAlphaHi*DToR));
CMCoef =(XM-0.25)/tan((VAlphaHi-90)*DToR);

%Create new table that cover entire AoA range.Fill in CL, CD, CM in original range
Count=1;
Step = 10;    %'angle of attack step in deg
SizeMax = 200;  %'dimension of Alpha arrays
HiEnd = 0;   %'intialize.  We haven't reached the hi end of the table
Remainder = VAlphaHi-floor(VAlphaHi/Step)*Step;

Alpha2(1) = -180;

for i=2:SizeMax
    Alpha2(i)= Alpha2(i-1)+Step;
    if Alpha2(i)>180-Step
        Alpha2(i) = 180;
        nTable2 = i;
        break
    end
    if Alpha2(i) >= VAlphaLo
        if Alpha2(i)-Step < VAlphaHi %use original data from 3DStall
            Alpha2(i) = Alpha1(Count);   
            CL2(i) = CL1(Count);
            CD2(i) = CD1(Count);
            CM2(i) = CM1(Count);
            Count=Count+1;
        elseif HiEnd == 0
            HiEnd = 1;
            Alpha2(i) = VAlphaHi - Remainder + Step;
        end
    end
end

%Function ViternaFill
CLAdj = 0.7;  %adjustment factor for CL when abs(Alpha)>90
for i=1:nTable2
    Alfa=Alpha2(i);
    SAlfa=sin(Alfa * DToR);
    CAlfa=cos(Alfa * DToR);
    
    if Alfa>180 || Alfa<-180
        warning('outside range + to -180 deg in Viterna calculation');
        break;
        
    elseif Alfa>=VAlphaHi && Alfa<= 90
        CL2(i)=CDMax/2*sin(2*Alfa*DToR)+ A2*CAlfa^2/SAlfa;
        CD2(i)=CDMax*SAlfa^2+B2*CAlfa;
        
    elseif Alfa>90 && Alfa<= 180-VAlphaHi     
        Ang=180-Alfa;
        Sang=sin(Ang*DToR);
        Cang=cos(Ang*DToR);
        CL2(i)=CLAdj*(-CDMax/2*sin(2*Ang*DToR)-A2*Cang^2/Sang);
        CD2(i)=CDMax*Sang^2+B2*Cang;
        
    elseif Alfa>180- VAlphaHi && Alfa<=180
        Ang=Alfa-180;
        Sang=sin(Ang*DToR);
        Cang=cos(Ang*DToR);
        CL2(i)=Ang/VAlphaHi*CLHi*CLAdj;
        CD2(i)=CDMax*Sang^2+B2*Cang;
        
    elseif Alfa>-VAlphaHi && Alfa< VAlphaLo
        CL2(i)=-CLHi*CLAdj+(Alfa+VAlphaHi)/(VAlphaHi+VAlphaLo)*(CLHi*CLAdj+CLLo); %Change by Andrew Ning.  The old code was: CL2(i) = CLAdj * (-CLHi + (Alfa + VAlphaHi) / (VAlphaHi + VAlphaLo) * (CLHi + CLLo))
        CD2(i)=CDLo+(-Alfa+VAlphaLo)*(CDHi-CDLo)/(VAlphaHi+VAlphaLo);

    elseif Alfa<=-VAlphaHi && Alfa>=-90
        Ang=-Alfa;
        Sang=sin(Ang*DToR);
        Cang=cos(Ang*DToR);
        CL2(i)=CLAdj*(-CDMax/2*sin(2*Ang*DToR)-A2*Cang^2/Sang);
        CD2(i)=CDMax*Sang^2+B2*Cang;
                
    elseif Alfa<-90 && Alfa>=-180+VAlphaHi
        Ang=180+Alfa;
        Sang=sin(Ang*DToR);
        Cang=cos(Ang*DToR);
        CL2(i)=CLAdj*(CDMax/2*sin(2*Ang*DToR)+A2*Cang^2/Sang);
        CD2(i)=CDMax*Sang^2+B2*Cang;

    elseif Alfa<-180+VAlphaHi && Alfa>=-180
        Ang=180+Alfa;
        Sang=sin(Ang*DToR);
        Cang=cos(Ang*DToR);
        CL2(i)=Ang/VAlphaHi*CLHi*CLAdj;
        CD2(i)=CDMax*Sang^2+B2*Cang;
    end
    
    if CD2(i)<0    %'watch out for negative CD's
        CD2(i)=0.01;
    end
end

%Do CM calculations and write value to output table
%Function GetCM
% x=coefficient used in CM calc.

 for i=1:nTable2
     Alfa = Alpha2(i);
     if Alfa>=VAlphaLo && Alfa<=VAlphaHi
        continue; %no action needed
     end
    
     if Alfa>-165 && Alfa<165
         if abs(Alfa)<0.01
             CM2(i)=CM0;
             
         elseif Alfa>0 
             x=CMCoef*tan((Alfa-90)*DToR)+0.25;
            CM2(i)=CM0-x*(CL2(i)*cos(Alfa*DToR)+CD2(i)*sin(Alfa*DToR));
         
         else 
             x=CMCoef*tan((-Alfa-90)*DToR)+0.25;
            CM2(i)=-(CM0-x*(-CL2(i)*cos(-Alfa*DToR)+CD2(i)*sin(-Alfa*DToR)));
         end
         
     else
         switch(Alfa)
             case 165
                 CM2(i)=-0.4;
             case 170
                 CM2(i)=-0.5;
             case 175
                 CM2(i) = -0.25;
             case 180
                 CM2(i) = 0;
             case -165
                 CM2(i) = 0.35;
             case -170
                 CM2(i) = 0.4;
             case -175
                 CM2(i) = 0.2
             case -180
                 CM2(i) = 0;
             otherwise
                 warning('Angle encountered for which there is no CM table value (near +/-180?).Program will stop');
         end     
     end  
 end  
Alpha2=Alpha2';
 CL2=CL2';
 CD2=CD2';
 CM2=CM2';
 CM0;
 nTable2;
 
end

