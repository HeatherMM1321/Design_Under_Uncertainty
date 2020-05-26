% Design Under Uncertainty : HW4 Utility Theory 
% Heather Miller 
% Started 5/25/20
% Due 5/27/20

clear all 
clc 

%adding CEtools folder to paths 
o_path = path;
path(o_path, 'C:\Users\MTH\Documents\Gradwork\Spring2020\Design_Under_Uncertainty\HW\HW4\CEtools');

problem = 4;

% used for problems 2-4 remaining problems 
samples = 10000;
rho = inf; %10 for risk adverse inf for risk neutral
high = 40;
low = 2; 

% Problem 1 - Determining the utility of A and D 
% uA~uB1+uC and uA~uB2+uD
if problem == 1
    prob_b_1 = 0.9;
    prob_b_2 = 0.7;
    prob_c = 0.1;
    prob_d = 0.3;

    utility_b = 1;
    utility_c = 0;

    utility_a = prob_b_1*utility_b + prob_c*utility_c;
    utility_d = (utility_a-(prob_b_2*utility_b))/prob_d;
    
    fprintf('Problem 1\n')
    fprintf('Utility of A: %.4f\n', utility_a)
    fprintf('Utility of D: %.4f\n', utility_d)
end 

% problem 2 - Monte Carlo method to compare the expected utility of 
% two different motor models 
if problem == 2
    d_1 = NaN(1, samples);
    d_2 = NaN(1, samples);
    ut_d1 = NaN(1, samples);
    ut_d2 = NaN(1, samples);
    
    for i = 1:samples
        %design variables 
        D = normrnd(7.5, 0.5);
        L = normrnd(9.5, 0.5);
        dcu = normrnd(8.94e3, 100);
        dfe = normrnd(7.98e3, 100);
        
        %weights given randomly generated values 
        d_1(i) = motorDesignHW4(D, L, dcu, dfe, 1);
        d_2(i) = motorDesignHW4(D, L, dcu, dfe, 2);
        
        %utility 
        ut_d1(i) = utility(d_1(i), rho, high, low);
        ut_d2(i) = utility(d_2(i), rho, high, low);
    end
    
    fprintf('Problem 2\n')
    fprintf('Design 1: %.4f\n', sum(d_1)/samples)
    fprintf('Utility 1: %.4f\n', sum(ut_d1)/samples)
    fprintf('Design 2: %.4f\n', sum(d_2)/samples)
    fprintf('Utility 2: %.4f\n', sum(ut_d2)/samples)
end

% problem 3 - Using FFNI method with 3 nodes to compare the expected 
% utility of 2 different motor models 
if problem ==3 
    indx = fullfact([3 3 3 3]);
    %design variables 
    [D, d_weight] = qnwnorm(3, 6.5, 0.5^2);
    [L, l_weight] = qnwnorm(3, 11.5, 0.5^2);
    [dcu, dcu_weight] = qnwnorm(3, 8.94e3, 100^2);
    [dfe, dfe_weight] = qnwnorm(3, 7.98e3, 100^2);
    
    design_1 = NaN(1, length(indx));
    design_2 = NaN(1, length(indx));
    W = NaN(1, length(indx));
    ut_d1 = NaN(1, length(indx));
    ut_d2 = NaN(1, length(indx));
    
    d1_CE = NaN(1, length(indx));
    d2_CE = NaN(1, length(indx));
    
    for i = 1: length(indx)
        % the index for each variable to use on this interation 
        D_idx = indx(i, 1);
        L_idx = indx(i, 2);
        dcu_idx = indx(i, 3); 
        dfe_idx = indx(i, 4); 
        %designs 
        design_1(i) = motorDesignHW4(D(D_idx), L(L_idx),dcu(dcu_idx),...
            dfe(dfe_idx), 1);
        design_2(i) = motorDesignHW4(D(D_idx), L(L_idx),dcu(dcu_idx),...
            dfe(dfe_idx), 2);
        % weights 
        W(i) = d_weight(D_idx)*l_weight(L_idx)*dcu_weight(dcu_idx)*....
            dfe_weight(dfe_idx);
        %utility 
        ut_d1(i) =  utility(design_1(i), rho, high, low)*W(i);
        ut_d2(i) =  utility(design_2(i), rho, high, low)*W(i);
        
        %CE design 1 
        d1_CE(i) = exp(design_1(i)/rho);
        
         %CE design 2 
        d2_CE(i) = exp(design_2(i)/rho);
        
       
    end
    % Design 1 moments 
    [d1_mean, d1_sigma, d1_skew, d1_kurt] = multi_moments(design_1, W);
    % Design 2 moments
    [d2_mean, d2_sigma, d2_skew, d2_kurt] = multi_moments(design_2, W);
    
    %CE design 1
    if rho == inf 
        d1_CE = CE_exponential(sum(design_1)/length(indx), rho);
    else 
        d1_CE = CE_exponential(sum(d1_CE)/length(indx), rho);
    end 
    
    if rho == inf 
        d2_CE = CE_exponential(sum(design_2)/length(indx), rho);
    else 
        d2_CE = CE_exponential(sum(d2_CE)/length(indx), rho);
    end 
   
    
    
    fprintf('Problem 3 : rho = %i\n', rho)
    fprintf('Design 1: %.4f\n', d1_mean)
    fprintf('CE 1: %.4f\n', d1_CE)
    fprintf('Utility 1: %.4f\n', sum(ut_d1))
    
    fprintf('Design 2: %.4f\n', d2_mean)
    fprintf('CE 2: %.4f\n', d2_CE)
    fprintf('Utility 2: %.4f\n', sum(ut_d2))
    
    
end 

% problem 4 - using Robust Design formulation and the FTNI with 5 nodes 
% to compare the certainty equivalent of two motor models 

if problem == 4 
    indx = fullfact([5 5 5 5]);
    
    %design variables 
    [D, d_weight] = qnwnorm(5, 6.5, 0.5^2);
    [L, l_weight] = qnwnorm(5, 11.5, 0.5^2);
    [dcu, dcu_weight] = qnwnorm(5, 8.94e3, 100^2);
    [dfe, dfe_weight] = qnwnorm(5, 7.98e3, 100^2);
    
    design_1 = NaN(1, length(indx));
    design_2 = NaN(1, length(indx));
    W = NaN(1, length(indx));
    ut_d1 = NaN(1, length(indx));
    ut_d2 = NaN(1, length(indx));
    
    for i = 1: length(indx)
        % the index for each variable to use on this interation 
        D_idx = indx(i, 1);
        L_idx = indx(i, 2);
        dcu_idx = indx(i, 3); 
        dfe_idx = indx(i, 4); 
        %designs 
        design_1(i) = motorDesignHW4(D(D_idx), L(L_idx),dcu(dcu_idx),...
            dfe(dfe_idx), 1);
        design_2(i) = motorDesignHW4(D(D_idx), L(L_idx),dcu(dcu_idx),...
            dfe(dfe_idx), 2);
        % weights 
        W(i) = d_weight(D_idx)*l_weight(L_idx)*dcu_weight(dcu_idx)*....
            dfe_weight(dfe_idx);
        %utility 
        ut_d1(i) =  utility(design_1(i), rho, high, low)*W(i);
        ut_d2(i) =  utility(design_2(i), rho, high, low)*W(i);
        
       
    end
    % Design 1 moments 
    [d1_mean, d1_sigma, d1_skew, d1_kurt] = multi_moments(design_1, W);
    % Design 2 moments
    [d2_mean, d2_sigma, d2_skew, d2_kurt] = multi_moments(design_2, W);
    
    %CE design 1
    d1_CE = d1_mean + ((d1_sigma^2)/(2*rho));
    d2_CE = d2_mean + ((d2_sigma^2)/(2*rho));
    
    
    fprintf('Problem 4: rho = %i\n', rho)
    fprintf('Design 1: %.4f\n', d1_mean)
    fprintf('CE 1: %.4f\n', d1_CE)
    fprintf('Utility 1: %.4f\n', sum(ut_d1))
    
    fprintf('Design 2: %.4f\n', d2_mean)
    fprintf('CE 2: %.4f\n', d2_CE)
    fprintf('Utility 2: %.4f\n', sum(ut_d2))
end 

