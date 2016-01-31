% SAVE HMM results in file /results/HMM.mat
% HMM transitions etc. are stored in the structure 
% HMMResults(i).change
% for each i=1,...,ntrials; 
%


HMM=1; % 1 to run HMM; 0 to load HMM
BICCRITERION=0;
AICCRITERION=1;
% 
NumFig=1;
close all;
addpath('../change_point');
addpath(genpath('../pmtk3-1nov12'));
% load data
fileload='../Data/LDdata/patientA.txt';
%fileload='Data/LDdata/patientB.txt';
Dataload=load(fileload);
seq=Dataload;

% figure folder
FigFolder='figs';
% output file
SaveFolder='results';


ntrials=10;
nspecies=size(seq,2);
%nspecies=300;

Lowpass=1;
% PARAMETERS
d=nspecies; % number of observations (emissions)
T = size(seq,1); % time steps

seq=seq(1:Lowpass:size(seq,1),:);

%id=find(seq==0);
sigma=min(seq(seq>0));
%seq(id)=5*rand(1,numel(id));
%seq=log(seq);
%
% REFORMAT seq
% Onset=5001;
% IN=input to CP function: structure array of length nspecies and field 'data'
% IN(t,n).data=array of dim: 1st=trials; 2nd=time bins
indSpecies=ntrials*ones(1,nspecies);
CumSpecies=[0 cumsum(indSpecies)];

 DATA=repmat(struct('data',[]),ntrials,nspecies);
app=repmat(struct('data',[]),ntrials,nspecies);
  for n=1:nspecies
    for t=1:ntrials
        %DATA(t,n).data=seq(:,CumSpecies(n)+t);
         DATA(t,n).data=seq(:,n)+0.01*sigma*rand(size(seq,1),1);
         app(t,n).data=seq(:,n);
    end
  end
 
% CREATE OBSERVATION SEQUENCE
observed=cell(ntrials,1);
for t=1:ntrials
    temp=zeros(nspecies,T);
    for n=1:nspecies
        temp(n,1:T)=app(t,n).data';
    end
    observed{t}=temp;
end

%----
% EM
%----
filename=fullfile(SaveFolder,'HMM_resultsLDA.mat'); %Dove salva i risultati
KSet=3;
if HMM
    Options.maxIter=1000; %Max Iteration
    Options.convTol=1e-7; %Diff tra due step EM
    Options.nRandomRestarts=7; %Expectation maximuzation start 10 times
    DATA=repmat(struct('modelEM',[],'loglikHist',[],'Posteriors',[]),1,numel(KSet));
%     parfor K_cnt=1:numel(KSet) % use this for parallelizing on numel(KSet) cores
    for K_cnt=1:numel(KSet)
        [DATA(K_cnt).modelEM, DATA(K_cnt).loglikHist, DATA(K_cnt).Posteriors]=HMMRun(observed,KSet(K_cnt),Options);
    end
    save(filename,'DATA');
else
    load(filename,'DATA');
end

    % BIC and AIC

    LL=zeros(1,numel(KSet));
    AIC=zeros(1,numel(KSet));
    BIC=zeros(1,numel(KSet));
    
    for K_cnt=1:numel(KSet)
        NParam=(d+d*(d+1)/2)*KSet(K_cnt)+ KSet(K_cnt)^2-KSet(K_cnt);
        LL(K_cnt)=-2*DATA(K_cnt).loglikHist(end);
        AIC(K_cnt)=-2*DATA(K_cnt).loglikHist(end)+2*NParam;
        BIC(K_cnt)=-2*DATA(K_cnt).loglikHist(end)+NParam*log(T);
    end

   % figure(NumFig); clf; NumFig=NumFig+1; hold on;
   % h(1)=plot(KSet,AIC,'b');
   % h(2)=plot(KSet,BIC,'r');
   % h(3)=plot(KSet,LL,'k');
   % legend(h,'AIC','BIC','LL');
   % saveas(gcf,fullfile('figs','LDB.BIC_VS_AIC.pdf'),'pdf');
   % hold off;
     if BICCRITERION
         [~,Crit]=min(BIC);
     elseif AICCRITERION
         [~,Crit]=min(AIC);
     end
 
    Detection=0.8;
    for Crit=1:numel(KSet)
        modelEM=DATA(Crit).modelEM;
        loglikHist=DATA(Crit).loglikHist;
        Posteriors=DATA(Crit).Posteriors;
        %ModelEM mi da concentrazione e varianza di ciascuna specie in ciascuno stato%
         % PARAMETERS
    
         %--------------------
         % EXTRACT TRANSITIONS
         %--------------------
        Options=[];
        Options.Detection=Detection;
        Options.BinSize=1;
        Sequences=HmmTransitions(Posteriors,Options); % tempo iniziale,tempo finale,durata,stato (righe) 
         %------------
        % PLOT TRIALS
         %------------
        % plot state sequence

        Options=[];
        Options.lowpass=1; % course grain 10 bins
        for i=1:1
             figure(NumFig); clf; NumFig=NumFig+1; hold on; 
             Posterior=Posteriors(i);
             Sequence=Sequences(i)
             Obs=observed{i};
             plotHMM(Obs,Posterior,Sequence,Options,Lowpass);
             filename=fullfile('figs',sprintf('LDA.HMM_trial%d_%d.pdf',KSet(Crit),i));
             saveas(gcf,filename,'pdf');
        end
    end

