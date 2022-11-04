% -------------------------------------------------------
%
%    f_getRRseries - This function generates an time series of
%    consecutive RR intervals
%
%    Ver. 1.0.0
%
%    Created:           Claudia Nagel (13.10.2022)
%    Last modified:     Claudia Nagel (13.10.2022)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2020 - All rights reserved.
%
% ------------------------------------------------------
%
% function [tau] = f_getRRseries(mu, N)
%
% Inputs:
%       mu: mean RR interval in ms
%       N: number or RR intervals to be generated
%
% Outputs:
%       tau: rr time series
%
%
% Example Usage:
%       rr = f_getRRseries(750, 10)
%
% Revision history:
%
%
function [tau] = f_getRRseries(mu, N)
% generates a time series of N beat-to-beat RR intervals based on a mean RR
% interval (mu) as proposed in Kantelhardt et al., 2002: "Modeling
% transient correlations in heartbeat dynamics during sleep"

% parameters taken from table I in Kantelhardt et al.:
gamma = 2.1;
b = 0.75; 

% initialize arrays
x = nan(N,1);
k = nan(N,1);
y = nan(N,1);
t = nan(N,1);
tau = nan(N,1);

% generate time series of N beats
for i = 1:N
    success = false; 
    while success == false
        x(i,1) = exprnd(1).*sign(rand-0.5); % Eq(1)

        k(i,1) = round((rand)^(1/(-gamma+1))); % Eq(2) & https://www.mathworks.com/matlabcentral/answers/294847-how-to-generate-power-law-random-numbers
        if i <= k(i,1)
            y(i,1) = 0; 
        else
            startI = i-k(i,1);
            endI = i-1;
            avrg_y = nanmean(y(startI:endI,1).^2);
            y(i,1) = x(i,1)*sqrt(1+b*avrg_y); % eq(3)
        end

        nbeats = 4;
        if i < nbeats+1
            t(i,1) = 0;
        else
            startI = i-mod(i,nbeats)-1;
            endI = i-1;
            t(i,1) = sum(tau(startI:endI,1));
        end
        yj_sum = 0; 
        for j = 1:i
            if j>i-k(j,1) && k(j,1)<i
                yj_sum = yj_sum + y(j,1); 
            end
        end
        tau(i,1) = mu+0.03*sin((2*pi*t(i,1))/3.6)+0.025*yj_sum;
        if tau(i,1)>2*mu || tau(i,1)<0.5*mu 
            success = false;
        else
            success = true; 
        end
    end
end
end



