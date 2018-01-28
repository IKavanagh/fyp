clear all
close all
clc


mu0 = 4*pi*1.0e-7 ; 
epsilon0 = 8.854*1.0e-12 ;
c0 = 1.0/sqrt(mu0*epsilon0) ; 

MHz = 1000000.0 ; 
frequency = 1000*MHz ; 
omega = 2.0*pi*frequency ;
lambda0 = c0/frequency ; 
k0 = omega/c0 ; 
eta1=sqrt(mu0/epsilon0);

small_threshold = 0.001 ; 

sigmac = 0.5  ;
muc = mu0 ;
epsilonc = 5*epsilon0 ; 

propagation_constant_c		=	sqrt( 1i*omega*muc*(sigmac +1i*omega*epsilonc) ) ;
kc = -1i*propagation_constant_c ; 
etac = 1i*omega*muc/propagation_constant_c ; 
lambdac = 2*pi/(real(kc)) ;
    
% epsilon = 5.0*epsilon0;
% c=1.0/sqrt(mu0*epsilon);e
k=kc ;
eta2= etac  ; 




n = 1.5;

xmin = -5 ; 
xmax = 0  ; 
ymin = -5 ; 
ymax = 0 ; 

% Building definition
X=[xmin xmax];
Y=[xmin xmax];
axis([-10 10 -10 10]);

x_edge1= xmin;
x_edge2= xmax;

h = 3 ; % Height of antenna above 0 face.  
TX=[xmin  h];  % TX location
hold on;
plot(TX(1,1),TX(1,2),'ko');

RX_x = 4.992515 ; 

h_dash = h*RX_x/(xmax - xmin) ;  % Used to compute where LOS and reflections occur. assumes vertical path. 


RX_Start=[RX_x -5]; %RX Line Start Point
RX_End=[RX_x  5];  %RX Line End Point
Line_Length = RX_Start(1,2)-RX_End(1,2);
Points = 2400;  % Total Points

step = Line_Length/Points;
RX_points = zeros(Points,2);

for (ct=1:Points+1)
    RX_points(ct,1) = RX_Start(1,1);
    RX_points(ct,2) = RX_Start(1,2)-((ct-1)*step);
    hold on
    plot(RX_points(ct,1),RX_points(ct,2),'k', RX_points(ct,1)+0.01,RX_points(ct,2),'k', RX_points(ct,1)-0.01,RX_points(ct,2),'k'); % Plot RX Line
end

s=zeros(Points,1);
LOS = zeros(Points,1) ; 

for(ct=1:Points)
    
    xs=RX_points(ct,1)-TX(1,1);
    ys=RX_points(ct,2)-TX(1,2);
    s(ct)=sqrt((xs*xs)+(ys*ys));    % Distance between TX and RX Point
    
    if( RX_points(ct,2) > -h_dash) 
        LOS(ct) = 1 ;
    
    hold on
    plot([TX(1,1) RX_points(ct,1)],[TX(1,2) RX_points(ct,2)],'r');  %Plot Direct Rays
    end
    
end

image = [-5 -h];            %Image Point
y_line = 0;     %visible edge line equation: y=3
ref_count=0;    %Reflection Points Counter 
reflection_point = zeros(Points,2);   %Reflection Points Vector
reflection_exists = zeros(Points,1) ; 

for (ct=1:Points)
    
    x_pt=(((y_line-RX_points(ct,2))*image(1,1))-((y_line-(image(1,2)))*RX_points(ct,1)))/(image(1,2)-RX_points(ct,2));  %See Report
    
    if( RX_points(ct,2) > h_dash )
        
    if(x_pt<=x_edge2 & x_pt>=x_edge1)
        reflection_exists(ct) = 1 ;  
        reflection_point(ct,1)=x_pt;
    end
    end
    
end

for (ct=1:ref_count)
    reflection_point(ct,2) = y_line;   %All reflection points lie on same edge
end

for(ct=1:Points)
hold on
if(reflection_exists(ct) == 1) 
plot([TX(1,1) reflection_point(ct,1) RX_points(ct,1)],[TX(1,2) reflection_point(ct,2) RX_points(ct,2)],'b');  
end

end

plot([X(1,1),X(1,2),X(1,2),X(1,1), X(1,1)],[Y(1,1),Y(1,1),Y(1,2),Y(1,2),Y(1,1)],'k');

s1=zeros(Points,1);
s2=zeros(Points,1);

for(ct=1:Points)

    xs1=reflection_point(ct,1)-TX(1,1);
    ys1=reflection_point(ct,2)-TX(1,2);
    s1(ct)=sqrt((xs1*xs1)+(ys1*ys1));   %Distances between TX and Reflection Point
    
    xs2=RX_points(ct,1)-reflection_point(ct,1);
    ys2=RX_points(ct,2)-reflection_point(ct,2);
    s2(ct)=sqrt((xs2*xs2)+(ys2*ys2));       %Distances between RX points and Reflection Point
end


theta_i=zeros(Points,1);
theta_t=zeros(Points,1);
Ref_Coefficient=zeros(Points,1);

for(ct=1:Points)%Calculation of Reflection Coefficients. See Report
    if(reflection_exists(ct) == 1) 
    perp=reflection_point(ct,1)-TX(1,1);
    theta_i(ct)=asin(perp/s1(ct));
    theta_t(ct)=asin((k0/k)*sin(theta_i(ct)));
    Ref_Coefficient(ct)=((eta2*cos(theta_i(ct)))-(eta1*cos(theta_t(ct))))/((eta2*cos(theta_i(ct)))+(eta1*cos(theta_t(ct))));
    end
end



    xs= x_edge2 - TX(1,1);
    ys= y_line - TX(1,2);
    dist_to_edge = sqrt(xs*xs + ys*ys) ;
    
    perp= x_edge2 - TX(1,1);
    
    theta_i_at_edge = asin(perp/dist_to_edge);
    theta_t_at_edge = asin((k0/k)*sin(theta_i_at_edge));
    Ref_Coefficient_at_edge = ((eta2*cos(theta_i_at_edge))-(eta1*cos(theta_t_at_edge)))/((eta2*cos(theta_i_at_edge))+(eta1*cos(theta_t_at_edge)));


E_inc=zeros(Points+1,1);    %Incident Field
E_ref=zeros(Points+1,1);    %Reflected Field
E_d1=zeros(Points+1,1);     %Diffracted field from Edge-1
E_d2=zeros(Points+1,1);     %Diffracted Field from Edge-2
E_GO=zeros(Points+1,1);     %Total GO Field (Inc + Ref)
E_D=zeros(Points+1,1);      %Total Diffracted(E_d1 + E_d2)
E_total=zeros(Points+1,1);  %Total Field (E_GO + E_D)
      
for(ct=1:Points)
  
   if(LOS(ct) == 1)
E_inc(ct)=besselh(0,2,k0*s(ct));    %Incident field in LOS points only
else
    E_inc(ct) = 0 ; 
end
    
end

for(ct=1:Points)
    if(reflection_exists(ct) == 1)
E_ref(ct)=(besselh(0,2,k0*s1(ct))*Ref_Coefficient(ct)*exp(-j*k0*s2(ct))*sqrt(s1(ct)/(s1(ct)+s2(ct))));   %Reflected Field

    end
end

    
for(ct=1:Points+1)
E_GO(ct)=E_inc(ct)+E_ref(ct);   %Total GO Field

end

%Diffraction from Edge-1. See Report

phi1_dash=pi-atan((TX(1,2)-y_line)/(x_edge1-TX(1,1)));
visible_pts_edge1=(Points/Line_Length)*(RX_Start(1,2)-y_line);  %Total visible points
rho1=zeros(visible_pts_edge1,1);
phi1=zeros(visible_pts_edge1,1);



 xs= x_edge1 - TX(1,1);
 ys= y_line - TX(1,2);   
dist_to_edge = sqrt(xs*xs + ys*ys) ; 
E_inc_at_edge = besselh(0,2,k0*dist_to_edge); 

for(ct=1:Points)
    
    if(RX_points(ct,2) > 0 ) 
        
    y21=RX_points(ct,2)-y_line;
    x21=RX_points(ct,1)-x_edge1;
   rho1(ct)= sqrt((y21*y21)+(x21*x21));
   phi1(ct)=atan(y21/x21);
   
   angle = phi1(ct) - phi1_dash ; 

   L = (rho1(ct) * dist_to_edge)/ ( rho1(ct) + dist_to_edge ) ; 

N_plus = round((pi + angle)/(2.0*n*pi)) ;
g_plus = 1.0 + cos(angle - 2.0*n*pi*N_plus) ;
%fresnel_plus = compute_fresnel(k0*rho1(ct)*g_plus) ; 
fresnel_plus = compute_fresnel(k0*L*g_plus) ; % Use L: Balanis p806.

N_minus = round((-pi + angle)/(2.0*n*pi)) ; 
g_minus = 1.0 + cos(angle - 2.0*n*pi*N_minus) ;

%fresnel_minus = compute_fresnel(k0*rho1(ct)*g_minus) ; 
fresnel_minus = compute_fresnel(k0*L*g_minus) ;  % Use L: Balanis p806.
   
D_i = -exp(-j*pi/4.0)/(2.0*n*sqrt(2.0*pi*k0))*( cot((pi + angle)/(2*n))*fresnel_plus + cot( (pi -angle)/(2*n))*fresnel_minus)  ;   
   
angle = phi1(ct) + phi1_dash ; 
N_plus = round((pi + angle)/(2.0*n*pi)) ;
g_plus = 1.0 + cos(angle - 2.0*n*pi*N_plus) ;
%fresnel_plus = compute_fresnel(k0*rho1(ct)*g_plus) ; 
fresnel_plus = compute_fresnel(k0*L*g_plus) ;  % Use L: Balanis p806.


N_minus = round((-pi + angle)/(2.0*n*pi)) ; 
g_minus = 1.0 + cos(angle - 2.0*n*pi*N_minus) ;
%fresnel_minus = compute_fresnel(k0*rho1(ct)*g_minus) ; 
fresnel_minus = compute_fresnel(k0*L*g_minus) ; 

D_r = -exp(-j*pi/4.0)/(2.0*n*sqrt(2.0*pi*k0))*( cot((pi + angle)/(2*n))*fresnel_plus+ cot( (pi -angle)/(2*n))*fresnel_minus )  ;
D_s = D_i  - D_r ;   %TMZ Polarization. Reference BALANIS 





E_d1(ct) = E_inc_at_edge*(exp(-j*k0*rho1(ct))/sqrt(rho1(ct)))*D_s ;   %Total Diffracted field from Edge-1
E_d1(ct) = 0 ; 

end
end

%Diffraction from Edge-2. See Report

phi2_dash = atan2((TX(1,2)-y_line),TX(1,1) - x_edge2 );
phi2_dash = pi - phi2_dash ; 


phi_2_dash = atan( (TX(1,2) - y_line)/(x_edge2 - TX(1,1) ) )  ; 


rho2=zeros(Points+1,1);
phi2=zeros(Points+1,1);

xs= x_edge2 - TX(1,1);
 ys= y_line - TX(1,2);


   
dist_to_edge = sqrt(xs*xs + ys*ys) ; 
E_inc_at_edge = besselh(0,2,k0*dist_to_edge);    

 
for(ct=1:Points)
    
%     if (ct>=visible_pts_edge1)      % For Points above Edge-2 height
%     y21=RX_points(ct,2)-y_line;
%     x21=RX_points(ct,1)-x_edge2;
%    rho2(ct)= sqrt((y21*y21)+(x21*x21));
%    phi2(ct)=pi-atan(y21/x21);
%     
%     elseif(ct==visible_pts_edge1+1)
%         rho2(ct)=RX_points(ct,1)-x_edge2;
%         phi2(ct)= pi;
%     else
%         y21=y_line-RX_points(ct,2);
%     x21=RX_points(ct,1)-x_edge2;
%    rho2(ct)= sqrt((y21*y21)+(x21*x21));
%    phi2(ct)=pi+atan(y21/x21);
%     end
%        

phi2(ct) =atan2((RX_points(ct,2)- y_line),(RX_points(ct,1)  - x_edge2 ));
phi2(ct) = pi - phi2(ct) ; 

phi2(ct) = pi + atan(( y_line - RX_points(ct,2)) / (RX_points(ct,1)  - x_edge2 ));


 y = RX_points(ct,2) - y_line;
 x= RX_points(ct,1) - x_edge2;
 rho2(ct)= sqrt((y*y)+(x*x));
   L = (rho2(ct) * dist_to_edge)/ ( rho2(ct) + dist_to_edge ) ; 
     
angle = phi2(ct) - phi2_dash ; 

N_plus = round((pi + angle)/(2.0*n*pi)) ;
g_plus = 1.0 + cos(angle - 2.0*n*pi*N_plus); 
%fresnel_plus = compute_fresnel(k0*rho2(ct)*g_plus);  
fresnel_plus = compute_fresnel(k0*L*g_plus);  

N_minus = round((-pi + angle)/(2.0*n*pi)) ;
g_minus = 1.0 + cos(angle - 2.0*n*pi*N_minus); 
%fresnel_minus = compute_fresnel(k0*rho2(ct)*g_minus) ;
fresnel_minus = compute_fresnel(k0*L*g_minus) ;
   
    D_i = -exp(-j*pi/4.0)/(2.0*n*sqrt(2.0*pi*k0))*( cot((pi + angle)/(2*n))*fresnel_plus + cot( (pi -angle)/(2*n))*fresnel_minus)  ;   
% else
% D_i = -1/(2.0)*(  sqrt(L)*sign(pi-angle) ) 
% end



angle = phi2(ct) + phi2_dash ; 
N_plus = round((pi + angle)/(2.0*n*pi)) ;
g_plus = 1.0 + cos(angle - 2.0*n*pi*N_plus) ;
%fresnel_plus = compute_fresnel(k0*rho2(ct)*g_plus) ; 
fresnel_plus = compute_fresnel(k0*L*g_plus) ; 

N_minus = round((-pi + angle)/(2.0*n*pi)) ; 
g_minus = 1.0 + cos(angle - 2.0*n*pi*N_minus) ;
%fresnel_minus = compute_fresnel(k0*rho2(ct)*g_minus) ; 
fresnel_minus = compute_fresnel(k0*L*g_minus) ; 

D_r = -exp(-j*pi/4.0)/(2.0*n*sqrt(2.0*pi*k0))*( cot((pi + angle)/(2*n))*fresnel_plus+ cot( (pi -angle)/(2*n))*fresnel_minus ) ;


D_s = D_i  + Ref_Coefficient_at_edge*D_r ;   %TMZ Polarisation 


E_d2(ct) = E_inc_at_edge*(exp(-j*k0*rho2(ct))/sqrt(rho2(ct)))*D_s ;   %Total Diffracted Field From Edge-2

end

for(ct=1:Points+1)
    E_D(ct)=E_d2(ct);          %Total Diffracted field(Edge-1 + Edge-2)
    E_total(ct)=E_GO(ct)+E_D(ct);       %Total Field(GO+Diffracted)
end


ct = 1:Points+1 ; 

figure
%plot(ct,(abs(E_inc(ct))),'r*');
%plot(ct,(abs(E_ref(ct))), 'g*');
plot(RX_points(ct,2),20*log10(abs(E_GO)),'r');
hold on
%plot(ct,(abs(E_d1(ct))),'r*');
%plot(ct,(abs(E_d2(ct))),'g*');
plot(RX_points(ct,2),20*log10(abs(E_D)),'g-.');
plot(RX_points(ct,2),20*log10(abs(E_total)),'k');

legend('E_G_O','E_D','Total');
xlabel(' y coord along vertical cut') 
ylabel(' E field (dB) ') 

 a= load('fields.res') ;
 
 
 figure
 plot(a(:,1), a(:,2), 'b-') ; 
 hold on
plot(RX_points(ct,2),20*log10(abs(E_total)),'r-.');
xlabel(' y coord along vertical cut') 
ylabel(' E  ') 
legend('VEFIE','UTD')  


