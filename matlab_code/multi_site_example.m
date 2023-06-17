cd 'F:\Work\RSPB\curlew_popModel\matlab_code'

% Matlab script for performing demographic, multi-site simulations.
% Morris & Doak 2002, Chapter 11, p435

%----- DemoMetaSim  -----
% Program to run multi-site stochastic simulations with correlations and
% a choice of different distributions for the vital rates.

clear all;

%-----  Simulation Parameters   -----
% order of the variables is survival and growth rates from the youngest to
% oldest, for each of 3 sites
% The fecundity used in all matrices:
%   s1a s2a s3a s4a g2a g3a s1b s2b s3b s4b
%   g2b g3b s1c s2c s3c s4c g2c g3c ff

load coryrates.txt;
% loads basic data with which to calculate all the relevant rates

vrates = coryrates;
% assign generic rate variable to these rates

vrates

% vrtypes identifies the distribution for each rate:
%   1 = beta
%   2 = stretched beta
%   3 = lognormal
% the function "ones" creates an array of all 1's for the first 13
% elements, followed by a single 3
% the number from the list for rate minimum and maximum values for each
% vital rate
% zeros are place holders for rates that are not stretched betas
vrtypes = [ones(1, 13), 3];

vrmins = zeros(1, 19);
vrmaxs = zeros(1, 19);

% create initial distribution of stages across sites
% the ' transposes from row to column
n0 = [720, 26, 13, 24, 920, 6, 0, 23, 980, 12, 20, 31]';

% Need a separate file that has MATLAB code that takes a vector of vital
% rate values and creates the matrix elements with them
% The program assumes the function uses a vector "vrs" to
% build the matrix
% assign matrix definition file: makes a separate matrix for each
%   pop. Essentially a grand matrix that contains the three single site
%       matrices
makemx = 'makecorymx';

% cap the populations of juveniles and adults
Ncap = [100, 100, 100];

% extinction threshold for each group
ne = [20 20 200];


tmax = 100;     % number of years to simulate
np = 19;        % number of vital rates
dims = 12;      % dimensions of the pop matrix
runs = 500;     % how many trajectories to do
popnum = 3;     % the number of separate pops
dimpop = 4;     % the number of classes per pop

%-----

vrmeans = mean(vrates);     % means of vital rates
vrvars = diag(cov(vrates)); % variances of vital rates
corrmatrix = corrcoef(vrates);  % correlation matrix for rates
                                % corrcoef refers to correlation coefs
corrmatrix(:, 19) = 0;      % the last param has no variance, so 
corrmatrix(19, :) = 0;      % correlation has NANs. Specific to this model
corrmatrix(19, 19) = 1;

for pp = 1:popnum   % total each pop size w/o seeds
    popNstart(pp) = sum(n0((2 + dimpop*(pp - 1)):dimpop*popnum));
end;

randn('state', sum(100 * clock)); % seeds random number generator (randn)
randn('state', sum(150 * clock));
[W, D] = eig(corrmatrix);   % make M12 matrix for correlations
M12 = W*sqrt(abs(D))*W';

%----- Make set of beta or stretched beta values to choose from  -----
%   make s99 values for 1% increments of Fx for each parameter
yesno = ...
    input('type 0 if you need to calculate betas, or 1 if you have a stored set');
if yesno == 1
    betafile = ...
        input('type filename with stored betas; put in single quotes');
    load betafile parabetas
else                % make a set of values for each beta or stretched beta
    parabetas = zeros(99, np+np2);
    for ii = 1:np
        if vrtypes(ii) ~= 3
            for fx99 = 1:99
                if vrtypes(ii) == 1;
                    parabetas(fx99, ii) = ...
                        betaval(vrmeans(ii), sqrt(vrvars(ii)), fx99/100);
                end;
                if vrtypes(ii) == 2
                    parabetas(fx99, ii) = ...
                        stretchbetaval(vrmeans(ii), sqrt(vrvars(ii)), ...
                        vrmins(ii), vrmaxs(ii), fx99/100);
                end;
            end;    %fx99
        end;        % if vrtypes
    end;            % ii loop

    yesno = input('type 1 to store the betas, or 0 if not');

    if yesno == 1
        betafile = ...
            input('type filename to store the betas, put in single quotes');
        save betafile, parabetas;
    end;            % if yesno
end;                % else

% the next section is same as 

