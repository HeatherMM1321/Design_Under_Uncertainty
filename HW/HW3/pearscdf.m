function [p,type] = pearscdf(X,mu,sigma,skew,kurt)

k=1;
plotf=1;
output=1;
method='G.Q.';
% pearspdf
%   [p,type,coefs] = pearspdf(X,mu,sigma,skew,kurt)
%
%   Returns the probability distribution denisty of the pearsons distribution
%   with mean `mu`, standard deviation `sigma`, skewness `skew` and
%   kurtosis `kurt`, evaluated at the values in X.
%
%   Some combinations of moments are not valid for any random variable, and in
%   particular, the kurtosis must be greater than the square of the skewness
%   plus 1.  The kurtosis of the normal distribution is defined to be 3.
%
%   The seven distribution types in the Pearson system correspond to the
%   following distributions:
%
%      Type 0: Normal distribution
%      Type 1: Four-parameter beta
%      Type 2: Symmetric four-parameter beta
%      Type 3: Three-parameter gamma
%      Type 4: Not related to any standard distribution.  Density proportional
%              to (1+((x-a)/b)^2)^(-c) * exp(-d*arctan((x-a)/b)).
%      Type 5: Inverse gamma location-scale
%      Type 6: F location-scale
%      Type 7: Student's t location-scale
%
%   Examples
%
%   See also
%       pearspdf pearsrnd mean std skewness kurtosis
%


%   References:
%      [1] Johnson, N.L., S. Kotz, and N. Balakrishnan (1994) Continuous
%          Univariate Distributions, Volume 1,  Wiley-Interscience.
%      [2] Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%          Springer-Verlag.

otpt=size(output,2);
outClass = superiorfloat(mu,sigma,skew,kurt);
%%
if X(2)==inf
    cdist = 1;
    limstate=X(1);
elseif X(1)==-inf
    cdist = 2;
    limstate=X(2);
else 
    cdist = 3;
    limstate=X;
end;

X = (X-mu)./sigma;  % Z-score

if strcmp(method,'MCS')
    beta1 = 0;
    beta2 = 3;
    beta3 = sigma.^2;
else
    beta1 = skew.^2;
    beta2 = kurt;
    beta3 = sigma.^2;
end;

% Return NaN for illegal parameter values.
if (sigma < 0) || (beta2 <= beta1 + 1)
    p = NaN(sizeOut,outClass);
    type = NaN;
    coefs = NaN(1,3,outClass);
    return
end

%% Classify the distribution and find the roots of c0 + c1*x + c2*x^2
c0 = (4*beta2 - 3*beta1); % ./ (10*beta2 - 12*beta1 - 18);
c1 = skew .* (beta2 + 3); % ./ (10*beta2 - 12*beta1 - 18);
c2 = (2*beta2 - 3*beta1 - 6); % ./ (10*beta2 - 12*beta1 - 18);
if c1 == 0 % symmetric dist'ns
    if beta2 == 3
        type = 0;
        a1=0;
        a2=0;
    else
        if beta2 < 3
            type = 2;
        elseif beta2 > 3
            type = 7;
        end
        a1 = -sqrt(abs(c0./c2));
        a2 = -a1; % symmetric roots
    end
elseif c2 == 0 % kurt = 3 + 1.5*skew^2 
    type = 3;
    a1 = -c0 ./ c1; % single root
    a2=a1;
else
    kappa = c1.^2 ./ (4*c0.*c2);
    if kappa < 0
        type = 1;
    elseif kappa < 1-eps
        type = 4;
    elseif kappa <= 1+eps
        type = 5;
    else
        type = 6;
    end
    % Solve the quadratic for general roots a1 and a2 and sort by their real parts
    tmp = -(c1 + sign(c1).*sqrt(c1.^2 - 4*c0.*c2)) ./ 2;
    a1 = tmp ./ c2;
    a2 = c0 ./ tmp;
    if (real(a1) > real(a2)), tmp = a1; a1 = a2; a2 = tmp; end
end

denom = (10*beta2 - 12*beta1 - 18);
if abs(denom) > sqrt(realmin)
    c0 = c0 ./ denom;
    c1 = c1 ./ denom;
    c2 = c2 ./ denom;
    coefs = [c0 c1 c2];
else
    type = 1; % this should have happened already anyway
    % beta2 = 1.8 + 1.2*beta1, and c0, c1, and c2 -> Inf.  But a1 and a2 are
    % still finite.
    coefs = Inf(1,3,outClass);
end

if strcmp(method,'MCS')
    type =8;
end;

%% Generate standard (zero mean, unit variance) values
switch type
case 0
    % normal: standard support (-Inf,Inf)
%     m1 = zeros(outClass);
%     m2 = ones(outClass);
    m1 = 0;
    m2 = 1;
    p = normcdf(X(2),m1,m2)-normcdf(X(1),m1,m2);
    Inv1=norminv( p, 0,1 );
%     Inv1=norminv( normcdf(X(1),m1,m2), 0,1 );    
    Inv2=norminv( normcdf(X(2),m1,m2), 0,1 );
    
    if plotf==1
        subplot(1,otpt,k)
        hold on;
        xx=(-6:0.01:6);
        col=[0 0 1];
        plot(xx*sigma+mu,(1/sigma)*normpdf(xx,m1,m2),'LineWidth',2)
        switch cdist
            case 1
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                line([limstate, limstate],[0,(1/sigma)*normpdf(X(1),m1,m2)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                line([limstate, limstate],[0,(1/sigma)*normpdf(X(2),m1,m2)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,(1/sigma)*normpdf(X(1),m1,m2)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,(1/sigma)*normpdf(X(2),m1,m2)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
        end;
        fill(xx1*sigma+mu,(1/sigma)*normpdf(xx,m1,m2),col)
        alpha(.2)
        title(output(k),'FontWeight','bold','FontSize',14);
        p_r=round(p*1000)/1000;
        height=norminv(.95,mu,sigma);
        text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
            if cdist==3 
                text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                'FontSize',14,'HorizontalAlignment','right')
                text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                'FontSize',14,'HorizontalAlignment','left')
            elseif cdist==2
                text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                'FontSize',14,'HorizontalAlignment','left')      
            else
                text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                'FontSize',14,'HorizontalAlignment','right')
            end;
    end;
case 1
        % four-parameter beta: standard support (a1,a2)
        if abs(denom) > sqrt(realmin)
            m1 = (c1 + a1) ./ (c2 .* (a2 - a1));
            m2 = -(c1 + a2) ./ (c2 .* (a2 - a1));
        else
            % c1 and c2 -> Inf, but c1/c2 has finite limit
            m1 = c1 ./ (c2 .* (a2 - a1));
            m2 = -c1 ./ (c2 .* (a2 - a1));
        end
    %     r = a1 + (a2 - a1) .* betarnd(m1+1,m2+1,sizeOut);
          X = (X-a1)./(a2-a1);        % Transform to 0-1 interval
    %     lambda = -(a2-a1)*(m1+1)./(m1+m1+2)-a1;
    %     X = (X - lambda - a1)./(a2-a1);
    disp(m1+1)
    disp(m2+1)

        p = betacdf(X(2),m1+1,m2+1)-betacdf(X(1),m1+1,m2+1);
        Inv1=norminv( p, 0,1 );
%         Inv1=norminv( betacdf(X(1),m1+1,m2+1), 0,1 );
        Inv2=norminv( betacdf(X(2),m1+1,m2+1), 0,1 );
        
        if plotf==1
            subplot(1,otpt,k)
            hold on;
            minp=(-6-a1)/(a2-a1);
            maxp=(6-a1)/(a2-a1);
            int=(maxp-minp)/1200;
            xx=(minp:int:maxp);
            plot((xx*(a2-a1)+a1)*sigma+mu,(1/(a2-a1))*(1/sigma)*(betapdf(xx,m1+1,m2+1)),'-r','LineWidth',2)
            switch cdist
            case 1
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                line([limstate, limstate],[0,(1/(a2-a1))*(1/sigma)*(betapdf(X(1),m1+1,m2+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                line([limstate, limstate],[0,(1/(a2-a1))*(1/sigma)*(betapdf(X(2),m1+1,m2+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,(1/(a2-a1))*(1/sigma)*(betapdf(X(1),m1+1,m2+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,(1/(a2-a1))*(1/sigma)*(betapdf(X(2),m1+1,m2+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;

            col=[1 0 0];
            fill((xx1*(a2-a1)+a1)*sigma+mu,(1/(a2-a1))*(1/sigma)*(betapdf(xx,m1+1,m2+1)),col)
            alpha(.2)
            title(output(k),'FontWeight','bold','FontSize',14);
            p_r=round(p*1000)/1000;
            height=norminv(.95,mu,sigma);
            text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                    'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
                if cdist==3 
                    text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                    text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                    'FontSize',14,'HorizontalAlignment','left')
                elseif cdist==2
                    text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                    'FontSize',14,'HorizontalAlignment','left')      
                else
                    text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                end;
        end;

%     X = X*(a2-a1) + a1;         % Undo interval tranformation
%     r = r + (0 - a1 - (a2-a1).*(m1+1)./(m1+m2+2));
case 2
    % symmetric four-parameter beta: standard support (-a1,a1)
    m = (c1+a1)./(c2.*2*abs(a1));
    m1=m;
    m2=m;
    X = (X-a1)./(2*abs(a1));
%     r = a1 + 2*abs(a1) .* betapdf(X,m+1,m+1);

        p = betacdf(X(2),m+1,m+1)-betacdf(X(1),m+1,m+1);
        Inv1=norminv( p, 0,1 );
%         Inv1=norminv( betacdf(X(1),m+1,m+1), 0,1 );
        Inv2=norminv( betacdf(X(2),m+1,m+1), 0,1 );

%     X = a1 + 2*abs(a1).*X;
        if plotf==1
            subplot(1,otpt,k)
            hold on;
            minp=(-6-a1)/(2*abs(a1));
            maxp=(6-a1)/(2*abs(a1));
            int=(maxp-minp)/1200;
            xx=(minp:int:maxp);
            plot((xx*(2*abs(a1))+a1)*sigma+mu,(1/(2*abs(a1)))*(1/sigma)*(betapdf(xx,m+1,m+1)),'-m','LineWidth',2)
            switch cdist
            case 1
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                line([limstate, limstate],[0,(1/(2*abs(a1)))*(1/sigma)*(betapdf(X,m+1,m+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                line([limstate, limstate],[0,(1/(2*abs(a1)))*(1/sigma)*(betapdf(X,m+1,m+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,(1/(2*abs(a1)))*(1/sigma)*(betapdf(X(1),m+1,m+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,(1/(2*abs(a1)))*(1/sigma)*(betapdf(X(2),m+1,m+1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;
            col=[1 0 1];
            fill((xx1*(2*abs(a1))+a1)*sigma+mu,(1/(2*abs(a1)))*(1/sigma)*(betapdf(xx,m+1,m+1)),col)
            alpha(.2)
            title(output(k),'FontWeight','bold','FontSize',14);
            p_r=round(p*1000)/1000;
            height=norminv(.95,mu,sigma);
            text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                    'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
                if cdist==3 
                    text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                    text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                    'FontSize',14,'HorizontalAlignment','left')
                elseif cdist==2
                    text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                    'FontSize',14,'HorizontalAlignment','left')      
                else
                    text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                end;
        end;
case 3
    % three-parameter gamma: standard support (a1,Inf) or (-Inf,a1)
    m = (c0./c1 - c1) ./ c1;
    m1=m;
    m2=m;
    X = (X - a1)./c1;
%     r = c1 .* gampdf(X,m+1,1,sizeOut) + a1;

        p = gamcdf(X(2),m+1,1)-gamcdf(X(1),m+1,1);
        Inv1=norminv( p, 0,1 );
%         Inv1=norminv( gamcdf(X(1),m+1,1), 0,1 );
        Inv2=norminv( gamcdf(X(2),m+1,1), 0,1 );

%     X = c1 .* X + a1;
        if plotf==1       
            subplot(1,otpt,k)
            hold on;
            minp=(-6-a1)/c1;
            maxp=(6-a1)/c1;
            int=(maxp-minp)/1200;
            xx=(minp:int:maxp);
            col=[1 .5 0];
            plot((xx*c1+a1)*sigma+mu,(1/c1)*(1/sigma)*(gampdf(xx,m+1,1)),'LineWidth',2,'LineStyle','-','Color',col)
            switch cdist
            case 1
                xx1=xx;
               xx1(xx1<=X(1))=X(1);
                line([limstate, limstate],[0,(1/c1)*(1/sigma)*(gampdf(X,m+1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                line([limstate, limstate],[0,(1/c1)*(1/sigma)*(gampdf(X,m+1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,(1/c1)*(1/sigma)*(gampdf(X(1),m+1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,(1/c1)*(1/sigma)*(gampdf(X(2),m+1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;   
            fill((xx1*c1+a1)*sigma+mu,(1/c1)*(1/sigma)*(gampdf(xx,m+1,1)),col)
            alpha(.2)
            title(output(k),'FontWeight','bold','FontSize',14);
            p_r=round(p*1000)/1000;
            height=norminv(.95,mu,sigma);
            text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                    'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
                if cdist==3 
                    text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                    text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                    'FontSize',14,'HorizontalAlignment','left')
                elseif cdist==2
                    text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                    'FontSize',14,'HorizontalAlignment','left')      
                else
                    text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                end;
        end;
case 4
    % Pearson IV is not a transformation of a standard distribution: density
    % proportional to (1+((x-lambda)/a)^2)^(-m) * exp(-nu*arctan((x-lambda)/a)),
    % standard support (-Inf,Inf)
     X=X.*sigma+mu;
      r   = 6*(beta2-beta1-1)/(2*beta2-3*beta1-6);
      m   = 1+r/2;
      nu  = -r*(r-2)*skew/sqrt(16*(r-1)-beta1*(r-2)^2);
      a   = sqrt(beta3*(16*(r-1)-beta1*(r-2)^2))/4;
      lambda = mu - ((r-2)*skew*sigma)/4        ;          % gives zero mean
      m1=m;
      m2=nu;
%     X = (X - lambda)./a;
        switch cdist
            case 1
                p = 1-pearson4cdf(X(1),m,nu,a,lambda,mu,sigma);
            case 2
                p = pearson4cdf(X(2),m,nu,a,lambda,mu,sigma);
            case 3
                p = pearson4cdf(X(2),m,nu,a,lambda,mu,sigma)-pearson4cdf(X(1),m,nu,a,lambda,mu,sigma);
        end;
        Inv1=norminv( p, 0,1 );
%         Inv1=norminv( pearson4cdf(X(1),m,nu,a,lambda,mu,sigma), 0,1 );
        Inv2=norminv( pearson4cdf(X(2),m,nu,a,lambda,mu,sigma), 0,1 );

%     C = X.*a + lambda;
%     C = diff(C);
%     C= C(1);
%     p = p./(sum(p)*C);
    if plotf==1
        subplot(1,otpt,k)
        hold on;
        minp=mu-6*sigma;
        maxp=mu+6*sigma;
        int=(maxp-minp)/1200;
        xx=(minp:int:maxp);
        col=[0 1 0];
        plot(xx,fx(xx,m,nu,a,lambda),'LineWidth',2,'Color',col)
            switch cdist
            case 1
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                line([limstate, limstate],[0,fx(X(1),m,nu,a,lambda)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                line([limstate, limstate],[0,fx(X(2),m,nu,a,lambda)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,fx(X(1),m,nu,a,lambda)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,fx(X(2),m,nu,a,lambda)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;  
        fill(xx1,fx(xx,m,nu,a,lambda),col)
        alpha(.2)
        title(output(k),'FontWeight','bold','FontSize',14);
        p_r=round(p*1000)/1000;
        height=norminv(.95,mu,sigma);
        text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
            if cdist==3 
                text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                'FontSize',14,'HorizontalAlignment','right')
                text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                'FontSize',14,'HorizontalAlignment','left')
            elseif cdist==2
                text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                'FontSize',14,'HorizontalAlignment','left')      
            else
                text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                'FontSize',14,'HorizontalAlignment','right')
            end;       
    end;
case 5
    % inverse gamma location-scale: standard support (-C1,Inf) or
    % (-Inf,-C1)
    C1 = c1 ./ (2*c2);
%     r = -((c1 - C1) ./ c2) ./ gampdf(X,1./c2 - 1,1) - C1;
    X = -((c1-C1)./c2)./(X + C1);
        m1=c2;
        m2=0;
        p = gamcdf(X(2),1./c2 - 1,1)-gamcdf(X(1),1./c2 - 1,1);
        Inv1=norminv( p, 0,1 );
%         Inv1=norminv( gamcdf(X(1),1./c2 - 1,1), 0,1 );
        Inv2=norminv( gamcdf(X(2),1./c2 - 1,1), 0,1 );

%     X = -((c1-C1)./c2)./X-C1;
        if plotf==1           
            subplot(1,otpt,k)
            hold on;
            minp=(-6-a1)/c1;
            maxp=(6-a1)/c1;
            int=(maxp-minp)/1200;
            xx=(minp:int:maxp);
            col=[1 .5 0];
            plot((xx*c1+a1)*sigma+mu,(1/c1)*(1/sigma)*(gampdf(xx,1./c2 - 1,1)),'LineWidth',2,'LineStyle','-','Color',col)
            switch cdist
            case 1
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                line([limstate, limstate],[0,(1/c1)*(1/sigma)*(gampdf(X,1./c2 - 1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                line([limstate, limstate],[0,(1/c1)*(1/sigma)*(gampdf(X,1./c2 - 1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,(1/c1)*(1/sigma)*(gampdf(X(1),1./c2 - 1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,(1/c1)*(1/sigma)*(gampdf(X(2),1./c2 - 1,1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;    
            fill((xx1*c1+a1)*sigma+mu,(1/c1)*(1/sigma)*(gampdf(xx,1./c2 - 1,1)),col)
            alpha(.2)
            title(output(k),'FontWeight','bold','FontSize',14);
            p_r=round(p*1000)/1000;
            height=norminv(.95,mu,sigma);
            text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                    'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
                if cdist==3 
                    text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                    text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                    'FontSize',14,'HorizontalAlignment','left')
                elseif cdist==2
                    text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                    'FontSize',14,'HorizontalAlignment','left')      
                else
                    text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                end;
        end;
case 6
    % F location-scale: standard support (a2,Inf) or (-Inf,a1)
    m1 = (a1 + c1) ./ (c2.*(a2 - a1));
    m2 = -(a2 + c1) ./ (c2.*(a2 - a1));
    % a1 and a2 have the same sign, and they've been sorted so a1 < a2
    if a2 < 0
        nu1 = 2*(m2 + 1);
        nu2 = -2*(m1 + m2 + 1);
        X = (X-a2)./(a2-a1).*(nu2./nu1);
%         r = a2 + (a2 - a1) .* (nu1./nu2) .* fpdf(X,nu1,nu2);

            p = fcdf(X(2),nu1,nu2)-fcdf(X(1),nu1,nu2);
            Inv1=norminv( p, 0,1 );
%             Inv1=norminv( fcdf(X(1),nu1,nu2), 0,1 );
            Inv2=norminv( fcdf(X(2),nu1,nu2), 0,1 );

%         X = a2 + (a2-a1).*(nu1./nu2).*X
        if plotf==1     
            subplot(1,otpt,k)
            hold on;
            minp=(-6-a2)/(a2-a1)*(nu2./nu1);
            maxp=(6-a2)/(a2-a1)*(nu2./nu1);
            int=(maxp-minp)/1200;
            xx=(minp:int:maxp);
            col=[.2 .3 .8];
            plot((xx*(a2-a1)/(nu2./nu1)+a2)*sigma+mu,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(xx,nu1,nu2)),'LineWidth',2,'LineStyle','-','Color',col)
            switch cdist
            case 1
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                line([limstate, limstate],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(1),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                line([limstate, limstate],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(2),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(1),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(2),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;    
            fill((xx1*(a2-a1)/(nu2./nu1)+a2)*sigma+mu,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(xx,nu1,nu2)),col)
            alpha(.2)
            title(output(k),'FontWeight','bold','FontSize',14);
            p_r=round(p*1000)/1000;
            height=norminv(.95,mu,sigma);
            text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                    'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
                if cdist==3 
                    text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                    text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                    'FontSize',14,'HorizontalAlignment','left')
                elseif cdist==2
                    text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                    'FontSize',14,'HorizontalAlignment','left')      
                else
                    text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                end;
        end;
    else % 0 < a1
        nu1 = 2*(m1 + 1);
        nu2 = -2*(m1 + m2 + 1);
        X = (X-a1)./(a1-a2).*(nu2./nu1);
%         r = a1 + (a1 - a2) .* (nu1./nu2) .* fpdf(X,nu1,nu2);
 
            p = -fcdf(X(2),nu1,nu2)+fcdf(X(1),nu1,nu2);
            Inv1=norminv( p, 0,1 );
%             Inv1=norminv( fcdf(X(1),nu1,nu2), 0,1 );
            Inv2=norminv( fcdf(X(2),nu1,nu2), 0,1 );

%         X = a1 + (a1-a2).*(nu1./nu2).*X;
        if plotf==1 
            subplot(1,otpt,k)
            hold on;
            minp=(-6-a1)/(a1-a2)*(nu2./nu1);
            maxp=(6-a1)/(a1-a2)*(nu2./nu1);
            int=(maxp-minp)/1200;
            xx=(minp:int:maxp);
            col=[.2 .3 .8];
            plot((xx*(a1-a2)/(nu2./nu1)+a1)*sigma+mu,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(xx,nu1,nu2)),'LineWidth',2,'LineStyle','-','Color',col)
            switch cdist
            case 1
                xx1=xx;
                xx1(xx1>=X(1))=X(1);
                line([limstate, limstate],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(1),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1<=X(2))=X(2);
                line([limstate, limstate],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(2),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1>=X(1))=X(1);
                xx1(xx1<=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(1),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(X(2),nu1,nu2))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;    
            fill((xx1*(a1-a2)/(nu2./nu1)+a1)*sigma+mu,((nu2./nu1)/(a2-a1))*(1/sigma)*(fpdf(xx,nu1,nu2)),col)
            alpha(.2)
            title(output(k),'FontWeight','bold','FontSize',14);
            p_r=round(p*1000)/1000;
            height=norminv(.95,mu,sigma);
            text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                    'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
                if cdist==3 
                    text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                    text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                    'FontSize',14,'HorizontalAlignment','left')
                elseif cdist==2
                    text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                    'FontSize',14,'HorizontalAlignment','left')      
                else
                    text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                    'FontSize',14,'HorizontalAlignment','right')
                end;
        end;
    end

case 7
    % t location-scale: standard support (-Inf,Inf)
    nu = 1./c2 - 1;
    X = X./sqrt(c0./(1-c2));
    m1=nu;
    m2=0;
        p = tcdf(X(2),nu)-tcdf(X(1),nu);
        Inv1=norminv( p, 0,1 );
%         Inv1=norminv( tcdf(X(1),nu), 0,1 );
        Inv2=norminv( tcdf(X(2),nu), 0,1 );

%     p = sqrt(c0./(1-c2)).*tpdf(X,nu);
%     X = sqrt(c0./(1-c2)).*X;
    if plotf==1    
        subplot(1,otpt,k)
        hold on;
        minp=-6/sqrt(c0./(1-c2));
        maxp=6/sqrt(c0./(1-c2));
        int=(maxp-minp)/1200;
        xx=(minp:int:maxp);
        col=[.9 .7 .7];
        plot(xx*sqrt(c0./(1-c2))*sigma+mu,(1/sigma)*(tpdf(xx,nu)),'LineWidth',2,'LineStyle','-','Color',col)
            switch cdist
            case 1
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                 line([limstate, limstate],[0,(1/sigma)*(tpdf(X,nu))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 2
                xx1=xx;
                xx1(xx1>=X(2))=X(2);
                 line([limstate, limstate],[0,(1/sigma)*(tpdf(X,nu))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            case 3
                xx1=xx;
                xx1(xx1<=X(1))=X(1);
                xx1(xx1>=X(2))=X(2);
                line([limstate(1), limstate(1)],[0,(1/sigma)*(tpdf(X(1),nu))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                line([limstate(2), limstate(2)],[0,(1/sigma)*(tpdf(X(2),nu))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;    

        fill(xx1*sqrt(c0./(1-c2))*sigma+mu,(1/sigma)*(tpdf(xx,nu)),col)
        alpha(.2)
        title(output(k),'FontWeight','bold','FontSize',14);
        p_r=round(p*1000)/1000;
        height=norminv(.95,mu,sigma);
        text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' method],  [' P=' num2str(p_r)]},...
                'FontSize',14,'HorizontalAlignment','left','EdgeColor',col,'Color',col)
            if cdist==3 
                text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
                'FontSize',14,'HorizontalAlignment','right')
                text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
                'FontSize',14,'HorizontalAlignment','left')
            elseif cdist==2
                text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
                'FontSize',14,'HorizontalAlignment','left')      
            else
                text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
                'FontSize',14,'HorizontalAlignment','right')
            end;
    end;
    
case 8
    %Monte Carlo Simulation Histogram
    out=kurt;
    p=skew;
    m1=0;
    m2=0;
    if plotf==1 
        subplot(1,otpt,k)
        hold on;
        nbins=25;
        [N,Xx]=hist(out,nbins);
        binwidth=Xx(2)-Xx(1);
        bar(Xx, N./(binwidth*sum(N)), 1,'cyan');
        title(output(k),'FontWeight','bold','FontSize',14);
            switch cdist
                case 1
                    line([limstate, limstate],[0,normpdf(limstate,mu,sigma)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                case 2
                    line([limstate, limstate],[0,normpdf(limstate,mu,sigma)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                case 3
                    line([limstate(1), limstate(1)],[0,normpdf(limstate(1),mu,sigma)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                    line([limstate(2), limstate(2)],[0,normpdf(limstate(2),mu,sigma)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
            end;
        p_r=round(p*1000)/1000;
        height=norminv(.95,mu,sigma);
        text(height,1.1*normpdf(height,mu,sigma),{['\leftarrow' 'MCS'],  [' P=' num2str(p_r)]},...
                'FontSize',14,'HorizontalAlignment','left','EdgeColor','cyan')
       if cdist==3 
            text(limstate(1),normpdf(limstate(1),mu,sigma),[num2str(limstate(1)) ' \rightarrow'],...
            'FontSize',14,'HorizontalAlignment','right')
            text(limstate(2),normpdf(limstate(2),mu,sigma),['\leftarrow' num2str(limstate(2)) ],...
            'FontSize',14,'HorizontalAlignment','left')
        elseif cdist==2
            text(limstate,normpdf(limstate,mu,sigma),['\leftarrow' num2str(limstate) ],...
            'FontSize',14,'HorizontalAlignment','left')      
        else
            text(limstate,normpdf(limstate,mu,sigma),[num2str(limstate) ' \rightarrow'],...
            'FontSize',14,'HorizontalAlignment','right')
       end;
    end;
end;

% scale and shift
% X = X.*sigma + mu; % Undo z-score
end

function p = pearson4cdf(X,m,nu,a,lambda,mu,sigma)
% pearson4pdf
%   p = pearson4pdf(X,m,nu,a,lambda)
%
%   Returns the pearson type IV probability density function with
%   parameters m, nu, a and lambda at the values of X.
%
%   Example
%
%   See also
%       pearson4pdf betapdf normpdf
%       pearspdf pearsrnd
%
Xx = (X - lambda)/a;

    if Xx<-sqrt(3)
        p1=fx(X,m,nu,a,lambda) *...
             a/(2*m-1) * (1i - Xx) *...
             hypergeom([1,m+nu/2*1i],2*m,2/(1-1i*Xx));
         p=real(p1);
    elseif Xx>sqrt(3)
        p1=1-fx(-X,m,-nu,a,-lambda) *...
             a/(2*m-1) * (1i + Xx) *...
             hypergeom([1,m-nu/2*1i],2*m,2/(1+1i*Xx));
        p=real(p1);
    elseif Xx<0 && Xx>-sqrt(3)&& abs(nu)<(4-2*sqrt(3))*m
%         p1=(1-exp(-(nu+1i*2*m)*pi))^(-1) -...
%            (1i*a*fx(X,m,nu,a,lambda))/...
%            (1i*nu-2*m+2) * (1+Xx^2) *...
%            hypergeom([1,2-2*m],2-m+1i*nu/2,(1+1i*Xx)/2);
        p1=normcdf(X,mu,sigma);
         p=real(p1);
    elseif Xx<0 && Xx>-sqrt(3)&& abs(nu)>(4-2*sqrt(3))*m
        p1=(1-exp(-(nu+1i*2*m)*pi))^(-1) -...
            (1i*a*fx(X,m,nu,a,lambda))/(1i*nu-2*m+2) * (1+Xx^2) *...
            hypergeom([1,2-2*m],2-m+1i*nu/2,(1+1i*Xx)/2);
        p=real(p1);
    elseif Xx>0 && Xx<sqrt(3)&& abs(nu)<(4-2*sqrt(3))*m
%         p1=1 - (1-exp(-(-nu+1i*2*m)*pi))^(-1) +...
%            (1i*a*fx(-X,m,-nu,a,-lambda))/...
%            (1i*(-nu)-2*m+2) * (1+(-Xx)^2) *...
%            hypergeom([1,2-2*m],2-m-1i*nu/2,(1-1i*Xx)/2);
       p1=normcdf(X,mu,sigma);
       p=real(p1);
    else
        p1=1 - (1-exp(-(-nu+1i*2*m)*pi))^(-1) +...
           (1i*a*fx(-X,m,-nu,a,-lambda))/...
           (1i*(-nu)-2*m+2) * (1+(-Xx)^2) *...
           hypergeom([1,2-2*m],2-m-1i*nu/2,(1-1i*Xx)/2);
        p=real(p1);
    end;

end

function [pearspdf] = fx(X,m,nu,a,lambda)
%Compute the pdf of pearson4
Xx = (X - lambda)/a;
k = -0.5.*log(pi)-log(a)-gammalnc(m-0.5)+...
     2*(real(gammalnc(m+(nu/2)*1i)))-gammalnc(m);
pearspdf =exp(k-m.*log(1+Xx.^2)-nu.*atan(Xx));

end

function [f] = gammalnc(z)
% GAMMALOG  Natural Log of the Gamma function valid in the entire complex plane.
%           This routine uses an excellent Lanczos series approximation
%           for the complex ln(Gamma) function.
%
%usage: [f] = gammaln(z)
%             z may be complex and of any size.
%             Also  n! = prod(1:n) = exp(gammalog(n+1))
%        
%tested under version 5.3.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%see also:   mhelp lnGAMMA

%Paul Godfrey
%pgodfrey@conexant.com
%07-13-01


siz = size(z);
z=z(:);
zz=z;

f = 0.*z; % reserve space in advance

p=find(real(z)<0);
if ~isempty(p)
   z(p)=-z(p);
end

%Lanczos approximation for the complex plane
 
g=607/128; % best results when 4<=g<=5
 
c = [  0.99999999999999709182;
      57.156235665862923517;
     -59.597960355475491248;
      14.136097974741747174;
      -0.49191381609762019978;
        .33994649984811888699e-4;
        .46523628927048575665e-4;
       -.98374475304879564677e-4;
        .15808870322491248884e-3;
       -.21026444172410488319e-3;
        .21743961811521264320e-3;
       -.16431810653676389022e-3;
        .84418223983852743293e-4;
       -.26190838401581408670e-4;
        .36899182659531622704e-5];

s=0;
for k=size(c,1):-1:2
    s=s+c(k)./(z+(k-2));
end

zg=z+g-0.5;
s2pi= 0.9189385332046727417803297;

f=(s2pi + log(c(1)+s)) - zg + (z-0.5).*log(zg);

f(z==1 | z==2) = 0.0;

if ~isempty(p)
   lpi= 1.14472988584940017414342735 + 1i*pi;
   f(p)=lpi-log(zz(p))-f(p)-log(sin(pi*zz(p)));
end

p=find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
   f(p)=Inf;
end

f=reshape(f,siz);
end

% Copyright (c) 2011, Christopher Hoyle.
% Developed with the sponsorship of the Defense Advanced Research Projects Agency (DARPA).

% Permission is hereby granted, free of charge, to any person obtaining a copy of this data, 
% including any software or models in source or binary form, as well as any drawings, 
% specifications, and documentation (collectively "the Data"), 
% to deal in the Data without restriction, including without limitation the rights to 
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Data, 
% and to permit persons to whom the Data is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Data.

% THE DATA IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS, SPONSORS, DEVELOPERS, CONTRIBUTORS, 
% OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
% WITH THE DATA OR THE USE OR OTHER DEALINGS IN THE DATA.