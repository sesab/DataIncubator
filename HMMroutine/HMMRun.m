% function [modelEM, loglikHist, Posterior]=HMMRun(observed,K)
% 
% HMM fit of observed sequences using EM
% INPUT: observed= cell array of length=# trials, each element is a matrix of dim d x T
%                 where d=# of species, T=length of trial
%        K       = number of hidden states
%        Options.maxIter=max # of EM iterations (default=1000)
%        Options.convTol=tolerance (default=1e-7)
%        Options.nRandomRestarts=# different EM fits (default=5)

% OUPUT: modelEM = structure with estimates from EM algorithm, fields:
%                 .nstates=K
%                 .type='gauss'
%                 .pi=probability of initial state
%                 .A=transition probability matrix
%                 .emission.mu=mean for each state
%                 .emission.Sigma=covariance matrix between emissions in each state
%                 etc.
%       loglikHist=log-likelihood at each time step
%       Posteriors=structure array of dim=# of trials, for each trial it has
%                 .viterbiPath=most likely sequence of states
%                 .gamma=posterior probability of each hidden state
%                 .loglik=loglikelihood
%                 .alpha, .beta=forward and backward probabilities
%                 .localEvidence=soft evidence
%                 .maxmargF=most likely state from alpha
%                 .maxmarg=most likely state from gamma
%
% Luca Mazzucato January 2016

                
function [modelEM, loglikHist, Posteriors]=HMMRun(observed,K,varargin)

% default
maxIter=1000;
convTol=1e-7;
nRandomRestarts=5;
if ~isempty(varargin)
    Options=varargin{1};
    if any(strcmp(fieldnames(Options),'maxIter'))
        maxIter=Options.maxIter;
    end
    if any(strcmp(fieldnames(Options),'convTol'))
        convTol=Options.convTol;
    end
    if any(strcmp(fieldnames(Options),'nRandomRestarts'))
        nRandomRestarts=Options.nRandomRestarts;
    end
end
%----
% EM
%----
DATA=repmat(struct('modelEM',[],'loglikHist',[]),1,nRandomRestarts);
%parfor cnt=1:nRandomRestarts
for cnt=1:nRandomRestarts
    % INITIALIZATION OF PRIORS
    DProb=0.990+0.01*(2*rand(K,1)-1); % random diagonal entries in initial t.p.m. in [0.98,1]
    trguess=zeros(K);
    for ent=1:K
        trguess(ent,1:K)=((1-DProb(ent))/(K-1))*ones(1,K);          % initial guess for the estimated probability of transition from state i to state j(ix)
        trguess(ent,ent)=DProb(ent);      % the probability of transitions are 0.001
    end

    
    nstates=K;
    [DATA(cnt).modelEM, DATA(cnt).loglikHist] = hmmFit(observed, nstates, 'gauss', ...
        'maxIter', maxIter, 'verbose', true, 'convTol', convTol, 'nRandomRestarts', 1,'transPrior',trguess);
end
LL=zeros(1,nRandomRestarts);
for cnt=1:nRandomRestarts
    LL(cnt)=DATA(cnt).loglikHist(end);
end
[~,maxcnt]=max(LL);
modelEM=DATA(maxcnt).modelEM;
loglikHist=DATA(maxcnt).loglikHist;

% PROBABILITIES
% logp = log p(X | model)
% alpha(i, t) = p(S(t)=i | X(:, 1:t)    (filtered)
% beta(i,t) propto p(X(:, t+1:T) | S(t=i))
% gamma(i,t)  = p(S(t)=i | X(:, 1:T))   (smoothed) % posterior probability
%   given whole sequence
% B - soft evidence
NSeq=numel(observed);
Posteriors=repmat(struct('viterbiPath',[],'gamma',[],'loglik',[],'alpha',[],'beta',[],...
    'localEvidence',[],'maxmargF',[],'maxmarg',[]),1,NSeq);
Names={'fieldnames','viterbiPath','gamma','loglik','alpha','beta',...
        'localEvidence','maxmargF','maxmarg'};
for i=1:NSeq
    Obs=observed{i};
%     len=size(Obs,2);
    %--------
    % VITERBI
    %--------
    viterbiPath = hmmMap(modelEM, Obs);
%     dgm = hmmToDgm(modelEM, len);
%     viterbiPathDGM = dgmMap(dgm, 'localev', Obs);
%     assert(isequal(viterbiPath, viterbiPathDGM));
%     % Sequence of Most Likely States (Max Marginals)
    [gamma, loglik, alpha, beta, localEvidence]  = hmmInferNodes(modelEM, Obs);
    % gamma is the posterior probability
    maxmargF = maxidx(alpha); % filtered (forwards pass only)
    maxmarg = maxidx(gamma);  % smoothed (forwards backwards)
    temp=v2struct(Names);
    Posteriors(i)=temp;
end
