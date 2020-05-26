function output = utility(x, rho, high, low)   
    if isinf(rho)
        output = (high-x)/(high-low);
    else
        output = (exp(-(high-x)/rho)-1)/(exp(-(high-low)/rho)-1);
        
    end
   
end
    
        