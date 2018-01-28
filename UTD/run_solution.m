% Run VEFIE solution for building
clear all
close all
clc

% Load variables
variables

% Create shape for problem to be solved over
create_shape

% Create VEFIE elements
create_vefie_elements

% Run VEFIE solution
vefie_solution



 the_max = -10000 ; 
    
    
    for(ct1=1:N)
        for(ct2=1:N)
    if(20*log10(abs(E_total(ct1,ct2))) > the_max ) 
        the_max = 20*log10(abs(E_total(ct1,ct2))) ; 
    end
        end
    end
    
 the_min = 10000 ; 
  
    for(ct1=1:N)
        for(ct2=1:N)
    if(20*log10(abs(E_total(ct1,ct2))) < the_min ) 
        the_min = 20*log10(abs(E_total(ct1,ct2))) ; 
    end
        end
    end
    
    for(ct1=1:N)
        for(ct2=1:N)
    
            if( building(ct1,ct2) == 1 ) 
         building(ct1,ct2) = the_max ; 
            end
            if( building(ct1,ct2) == 0 )
         building(ct1,ct2) = the_min ; 
    end
    
        end
    end
    
    
    
% Plot results
    h1 = figure;
    
    surf(real(pos), imag(pos), 20.0*log10(abs(E_total)));
    
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Electric Field throughout building');

    h1 = figure;
    
    surf(real(pos), imag(pos), log10(abs(E_inc)));
    
    hold on
    
    surf(real(pos), imag(pos), building);
    
    shading interp;
    view(2);
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    title('Incident Electric Field throughout building');

    
    figure
    plot( imag(pos(1:N,N) ) , 20*log10(abs(E_total(1:N,N) ) ))
    %legend('horizontal','vertical') ;
    xlabel(' y coord ' ) 
    ylabel(' field (dB) ')  ;
    
    fp = fopen('fields.res','w') ;
    
    for( ct = 1:N) 
    fprintf(fp,' %f %f \n',imag(pos(ct,N)),20*log10(abs(E_total(ct,N))) ) ; 
    
    end
    
    fclose(fp) ; 
    
    
    fprintf(1,'Computing fields along line x = %f \n',real(pos(1,N)) ) ; 
    
    
