% function CircIndent2LoadProfile
%
% compute load profile from displacement as in SkinTouchProfile.m
%
% Inputs : S0, indentation matrix (nsamples,npins)
%          xy, 2D coordinates of each pins (npins,2)
%
% Ouputs : P, equivalent load profile at the same coordinates
%          Pdyn, dynamics load profile
%          S1, corrected indentation matrix
%
% Notes: E=50kPa; nu=0.4;

function [P_Old, Pdyn_Old, S1_old] = CircIndent2LoadProfile_Old(S0,xy,samp_freq,ProbeRad)

s=size(S0);

if(s(2)~=size(xy,1))
    disp(['Trace size : ' num2str(size(S0))])
    disp(['Pin number : ' num2str(size(xy,1))])
    error('incoherent number of pins ');
end

% Young's Modulus
E=0.035000; % 50kPa = 50,000 N/m^2 = 0.05 N/mm^2;
nu=0.4;

x=xy(:,1);
y=xy(:,2);
% distance matrix
R=sqrt((bsxfun(@minus,x,x')).^2 + (bsxfun(@minus,y,y')).^2);

% flat cylinder indenter solution from (SNEDDON 1946)
D = (1-nu^2)/pi/ProbeRad*asin(ProbeRad./R)/E;
D(R<=ProbeRad)=(1-nu^2)/2/ProbeRad/E;

% hack to specify stretch (for glued probe)
% DOES NOT WORK IF A STIMULUS HAS POSITIVE AND NEGATIVE VALUES IN S0
% NON CONTACTING PINS SHOULD BE SET TO 0
% I.E. YOU CAN ONLY SPECIFY SKIN STRESHING/PULLING (S0<0) FOR FLAT SURFACE
% INDENTATIONS
if(any(S0<0))
    warning('Negative indentation found: proceed with caution !')
end

S0neg=S0<0;
absS0=abs(S0);

P_Old = zeros(s); prevS0=zeros(s); count=0;
% iterative contact-detection algorithm
while(count==0||~isempty(find(P_Old<0,1)))
    absS0(P_Old<0) = 0;
    count=count+1;
    % only work on changed (and nonzeros) line
    diffl=sum(absS0-prevS0,2)~=0;
    S0loc=absS0(diffl,:);
    P_Old(diffl,:)=blockSolve(S0loc,D); % efficiently solves P=D\S0loc ignoring zero-entry line in S0loc
    prevS0=absS0;
end

% correct for the hack
P_Old(S0neg)=-P_Old(S0neg);

% actual skin profile under the pins
if(nargout>1)
    S1_old=P_Old*D;
end

% time derivative of deflection profile
% assumes same distritubution of pressure as in static case
% proposed by BYCROFT (1955) and confirmed by SCHMIDT (1981)
opts.SYM = true; opts.POSDEF=true;
if(nargout>1)
    if(s(1)>1)
        % compute time derivative
        S1p=([S1_old(2:end,:); nan(1,size(S1_old,2))] - [nan(1,size(S1_old,2)) ; S1_old(1:end-1,:)])/2*samp_freq;
        S1p(1,:)=S1p(2,:); S1p(end,:)=S1p(end-1,:);
        % linsolve
        Pdyn_Old=linsolve(D,S1p',opts)'/1;
    else
        Pdyn_Old=zeros(size(P_Old));
    end
end
end

% this blockSolve function computes the operation P=D\S0 for non-zeros lines
% and for every time sample but avoids to loop through all time samples by
% indentifying similar blocks (blocks having the same non-zero patterns
% in S0 across time) and solves the equation for each block (instead
% of each line).
function P_Old=blockSolve(S0,D)
nz=S0~=0; % non zeros elements
% find similar lines to solve the linear system
B=bwpack(nz')';
[~,ia,ic]=unique(B,'rows'); % unique is much faster after bwpack (I guess booleans are uint32 or so)
unz=nz(ia,:); % unique non-zeros elements
opts.SYM = true; opts.POSDEF=true;
P_Old=zeros(size(S0));

for ii=1:length(ia)
    lines=(ic==ii);    % lines of this block
    nzi=unz(ii,:);     % non-zeros elements
    P_Old(lines,nzi) = linsolve(D(nzi,nzi),S0(lines,nzi)',opts)';
end
end