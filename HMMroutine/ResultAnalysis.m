% SAVE HMM results in file /results/HMM.mat
% HMM transitions etc. are stored in the structure 
% HMMResults(i).change
% for each i=1,...,ntrials; 
%
% Luca Mazzucato December 2016

HMM=1; % 1 to run HMM; 0 to load HMM
BICCRITERION=0;
AICCRITERION=1;
% 
NumFig=1;
close all;
%addpath('..change_point');
addpath(genpath('../pmtk3-1nov12'));
% load data
%fileload='../Data/LDdata/patientA1.txt';
fileload='../Data/LDdata/patientB.txt';
%fileload='../Data/Mouse/RDMouse.txt';
file1='PatientB';
SaveFolder='results';
FigFolder='figs';
KSet=3;
filename=fullfile(SaveFolder,sprintf('%s.HMM_result%d.mat',file1,KSet)) %Dove salva i risultati

Dataload=load(fileload);
seq=Dataload;

% PARAMETERS


ntrials=1;
nspecies=size(seq,2);
%nspecies=13;
Lowpass=1;

d=nspecies; % number of observations (emissions)
T = size(seq,1); % time steps

seq=seq(1:Lowpass:size(seq,1),:);

%id=find(seq==0);
sigma=min(seq(seq>0));
%Remember to remove 0
%seq(seq==0)=sigma*.01;
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
        
        %THIS IS FOR DATA THAT HAS MORE THAN ONE TRIALS
        %DATA(t,n).data=seq(:,CumSpecies(n)+t);
        %app(t,n).data=seq(:,CumSpecies(n)+t);
        %app(t,n).data=seq(:,CumSpecies(n)+t);
        
        temp=seq(:,n)+0.001*sigma*rand(size(seq,1),1);
      %  if(temp<0)
        % DATA(t,n).data = temp;
       % else
         DATA(t,n).data = temp;
       % end
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


if HMM
    Options.maxIter=1000; %Max Iteration
    Options.convTol=1e-5; %Diff tra due step EM
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

    %PLOT LL, AIC and BIC
    %figure(NumFig); clf; NumFig=NumFig+1; hold on;
    %h(1)=plot(KSet,AIC,'b');
    %h(2)=plot(KSet,BIC,'r');
    %h(3)=plot(KSet,LL,'k');
    %legend(h,'AIC','BIC','LL');
    %saveas(gcf,fullfile('figs','LDB.BIC_VS_AIC.pdf'),'pdf');
    %hold off;
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
        Options.MinDur=0.02;    %default .01*T
        Sequences=HmmTransitions(Posteriors,Options); % tempo iniziale,tempo finale,durata,stato (righe) 
         %------------
        % PLOT TRIALS
         %------------
        % plot state sequence

       
    end
    
    
    cnt=1; % Select the first KSet(cnt)
    if(numel(KSet)>1)
       fprintf('ResultAnalysis: %d possible selection of hidden states. We pick the smallest\n',KSet);
    end
    Cov1=DATA(cnt).modelEM.emission.Sigma(:,:,1);
    sigma=diag(1./sqrt(diag(Cov1)));
    Corr=sigma*Cov1*sigma;
    figure(10)
    J=inv(Corr);
    %[J1,id]=sort(J(:),'descend');
    
   imagesc(J)
   
      Cov2=DATA(cnt).modelEM.emission.Sigma(:,:,2);
    sigma=diag(1./sqrt(diag(Cov2)));
    Corr=sigma*Cov2*sigma;
    figure(11)
    imagesc(inv(Corr))
    
   if(KSet(cnt)>2)  
    Cov3=DATA(cnt).modelEM.emission.Sigma(:,:,3);
    sigma=diag(1./sqrt(diag(Cov3)));
    Corr=sigma*Cov3*sigma;
    figure(12)
    imagesc(inv(Corr))
   end
   
   
    
    
    Average1=DATA(cnt).modelEM.emission.mu(:,1);
    Average2=DATA(cnt).modelEM.emission.mu(:,2);
    if(KSet(cnt)>2)
      Average3=DATA(cnt).modelEM.emission.mu(:,3);
    end
    
    figure(13)
    y1=(Average2-Average1)./sqrt(diag(Cov1)+diag(Cov2))*sqrt(nspecies);
    if(KSet(cnt)>2)
     y2=(Average3-Average1)./sqrt(diag(Cov1)+diag(Cov3))*sqrt(nspecies);
    
     y3=(Average3-Average2)./sqrt(diag(Cov2)+diag(Cov3))*sqrt(nspecies);
    end
    [y1,id]=sort(y1,'descend');
   
    
    %%%%% COMPUTE CHI SQUARE FOR ALL POSSIBLE STATES
    L=numel(KSet)*(numel(KSet)-1)/2;
    stat=repmat(struct('h',[],'chi2Stat',[],'pvalue',[]),1,L);
    
    for i=1:KSet(cnt)
        mu1=DATA(cnt).modelEM.emission.mu(:,i);
        cov1=diag(DATA(cnt).modelEM.emission.Sigma(:,:,i));
        for j=1:KSet(cnt)
            mu2=DATA(cnt).modelEM.emission.mu(:,j);
            cov2=diag(DATA(cnt).modelEM.emission.Sigma(:,:,j));
            if(i<j)
               k=(i-1)*KSet(cnt)+j;
              
                [stat(k).h,stat(k).chi2Stat,stat(k).pvalue]=chi2Fun(mu1,cov1,mu2,cov2);
                distance(k)=sum(stat(k).h)/numel(stat(k).h);
            else
                k=(i-1)*KSet(cnt)+j;
                distance(k)=-1.;
            end
        end
    end
    
    
    if(KSet(cnt)==2)
      [dis,id]=sort(distance,'descend');
        k=id(1);
        [y,id1]=sort(stat(k).pvalue);
        
        bar(1:nspecies,-log(y),'b');
        hold on
        v=-log(.05)*ones(1,nspecies);
        plot(1:nspecies,v,'g','LineWidth',2);
        xlim([1 nspecies]);
        %ylim([min(y) 1]);
         [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
    end
    
    if(KSet(cnt)==3)
        [dis,id]=sort(distance,'descend');
        k=id(1);
        [y,id1]=sort(stat(k).pvalue);
        subplot(2,2,1)
        bar(1:nspecies,-log(y),'b');
        hold on
         v=-log(.05)*ones(1,nspecies);
        plot(1:nspecies,v,'g','LineWidth',2);
       
        xlim([1 nspecies]);
       % ylim([min(y) 1]);
        [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
        
    
        subplot(2,2,2)
        k=id(2);
        bar(1:nspecies,-log(stat(k).pvalue(id1)),'r');
        hold on
        xlim([1 nspecies]);
        
        plot(1:nspecies,v,'g','LineWidth',2);
        [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
    %legend('state %d and %d %lf',k/KSet(cnt),mod(k,K(set)),distance(k));
   
        k=id(3);
        subplot(2,2,[3 4])
        bar(1:nspecies,-log(stat(k).pvalue(id1)),'k');
        hold on
        
        plot(1:nspecies,v,'g','LineWidth',2);
        xlim([1 nspecies]);
       [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
    %legend('state 3-state 1')
    end
    
    
      if(KSet(cnt)==4)
        [dis,id]=sort(distance,'descend');
        k=id(1);
        [y,id1]=sort(stat(k).pvalue);
        subplot(2,3,1)
        bar(1:nspecies,-log(y),'b');
        hold on
         v=-log(.05/nspecies)*ones(1,nspecies);
        plot(1:nspecies,v,'g','LineWidth',2);
       
        xlim([1 nspecies]);
       % ylim([min(y) 1]);
        [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
        xlabel('Species')
    ylabel('-log Pvalue');
        
    
        subplot(2,3,2)
        k=id(2);
        bar(1:nspecies,-log(stat(k).pvalue(id1)),'r');
        hold on
        xlim([1 nspecies]);
        
        plot(1:nspecies,v,'g','LineWidth',2);
        [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
        xlabel('Species')
    ylabel('-log Pvalue');
    %legend('state %d and %d %lf',k/KSet(cnt),mod(k,K(set)),distance(k));
   
        k=id(3);
        subplot(2,3,3)
        bar(1:nspecies,-log(stat(k).pvalue(id1)),'k');
        hold on
        
        plot(1:nspecies,v,'g','LineWidth',2);
        xlim([1 nspecies]);
       [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
        xlabel('Species')
    ylabel('-log Pvalue');
    %legend('state 3-state 1')
        k=id(4);
        subplot(2,3,4)
        bar(1:nspecies,-log(stat(k).pvalue(id1)),'k');
        hold on
        
        plot(1:nspecies,v,'g','LineWidth',2);
        xlim([1 nspecies]);
       [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
        xlabel('Species')
    ylabel('-log Pvalue');
        
        k=id(5);
        subplot(2,3,5)
        bar(1:nspecies,-log(stat(k).pvalue(id1)),'c');
        hold on
        
        plot(1:nspecies,v,'g','LineWidth',2);
        xlim([1 nspecies]);
       [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
        xlabel('Species')
        ylabel('-log Pvalue');
   
         k=id(6);
        subplot(2,3,6)
        bar(1:nspecies,-log(stat(k).pvalue(id1)),'y');
        hold on
        
        plot(1:nspecies,v,'g','LineWidth',2);
        xlim([1 nspecies]);
        [i,j]=ind2sub(KSet,k);
        str=sprintf('%d and %d dist: %2.1f%%',i,j,distance(k)*100);
        legend(str,'Location','northwest');
        xlabel('Species')
    ylabel('-log Pvalue');
    end

    
    
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [50 0 700 1024]/100);
    filename=fullfile('figs',sprintf('%s.pvalue.species%d_%d.pdf',file1,KSet(cnt),i));
    print('-dpdf','-r200',filename); 
    
    
    Options=[];
    Options.lowpass=1; % course grain 10 bins
    Options.id=id1;  % just add for Mouse thick per species
    
     figure(NumFig); clf; NumFig=NumFig+1; hold on; 
    Posterior=Posteriors(cnt);
    Sequence=Sequences(cnt);
    i=1;
    Obs=observed{i};
     %plotHMMMouse(Obs,Posterior,Sequence,Options,Lowpass);
    plotHMM(Obs,Posterior,Sequence,Options,Lowpass);
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [50 0 700 1024]/100);

   filename=fullfile('figs',sprintf('%s.HMM_trial%d_%d.pdf',file1,KSet(Crit),i));
      print('-dpdf','-r200',filename); 
%              saveas(gcf,filename,'pdf');
        
   

   