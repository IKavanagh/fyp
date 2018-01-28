test_y_line = N/2  + 10 ; 
ct1 = test_y_line ; 

face_x = wall_start  - 0.5*delta_x ; 

wall_width = (w_dx+1)*delta_x ; 

R1 = (etac - eta0)/(eta0 + etac) ; 
R2 = (eta0 - etac)/(eta0 + etac); 
T1 = 1 + R1 ; 
T2 = 1 + R2 ;

a1 = exp(-2*1i*kc*(wall_width-face_x))*R2*R2 ; 
a2 = R2*R2*exp(-2*1i*kc*(wall_width-face_x)) ; 

a1 = exp(-2*1i*kc*(wall_width))*R2*R2 ; 
a2 = R2*R2*exp(-2*1i*kc*(wall_width)) ; 

analytic_E = zeros(N,1) ;
analytic_E_single = zeros(N,1) ;     

    field_at_face = exp(-1i*k0*face_x) ; 
    
for(ct3 = 1:N) 

    location(ct3) = start_point + ((ct3 - 0.5)*delta_x + (ct1 - 0.5)*delta_y*1i) ;
    the_x = real(location(ct3)) ; 
    the_y = imag(location(ct3)) ; 
    

    term1 = field_at_face*T1*T2*exp(-1i*k0*(the_x - (face_x + wall_width) ))*exp(-1i*kc*(wall_width)) ;
    term2 = field_at_face*exp(1i*k0*(the_x - face_x) )*exp(-1i*kc*2*(wall_width))*T1*T2*R2 ;
    
    
    
    if( the_x > wall_width) 
        
    analytic_E(ct3) = term1/(1-a1) ; 
    analytic_E_single(ct3) = term1 ; 
    end
    
    if( the_x < 0.0) 
    analytic_E(ct3) = exp(-1i*k0*the_x) + field_at_face*exp(1i*k0*(the_x - face_x))*R1 ; 
    analytic_E_single(ct3) = analytic_E(ct3) ;
    analytic_E(ct3) = analytic_E(ct3) +  term2/(1-a2) ; 
    end
    
end



figure 
plot(real(location),20*log10(abs(analytic_E_single))) 
hold on
plot(real(location),20*log10((analytic_E)),'k') 
%plot(real(location),real(E_inc(ct1,:)),'r')
plot(real(location),20*log10(abs(E_total(ct1,:))),'r')
xlabel('x coord')
ylabel('E') 
legend('single','ray','VEFIE')
        
  