function certainty_equivalence = CE_exponential(x, rho)
    if isinf(rho)
        certainty_equivalence = x;
    else
        certainty_equivalence = rho*log(x);
        
    end
   
end
