function result = compute_fresnel(argument) 
result = 0.0  ;
fresnel_data ; 

if(argument <= 0.3 )
result = (sqrt(pi*argument) - 2.0*argument*exp(j*pi/4.0) - (2/3)*power(argument,2)*exp(-j*pi/4.0))*exp(j*(pi/4.0 + argument)); 
end

if(argument >= 5.5)
result = (1.0 + j/(2.0*argument) - 0.75/(power(argument,2.0))  - (15*j/8)/power(argument,3) + (75/16)/(power(argument,4.0)) ) ; 
end

if( (argument > 0.3) && (argument < 5.5) ) 

index = floor(sqrt(argument)/0.1)   ; 

factor1 = sqrt(argument) - index*0.1 ; 
factor2 = 0.1 - factor1 ;

result = (C1(index+1)*factor2 + C1(index+2)*factor1)/0.1 ; 
result = result - j*(S1(index+1)*factor2 + S1(index+2)*factor1)/0.1 ; 
result = result*2.0*j*sqrt(argument)*exp(j*argument) ; 

end


