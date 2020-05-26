function [mean, sigma, skew, kurt] = multi_moments(output, Weight)
    mean = dot(Weight, output);
    sigma = sqrt(dot(Weight, (output-mean).^2));
    skew = dot(Weight, (output-mean).^3)/sigma.^3;
    kurt = dot(Weight, (output-mean).^4)/sigma.^4;
end
