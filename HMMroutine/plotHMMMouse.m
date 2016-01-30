function plotHMMMouse(y,Posterior,Sequence,Options,Lowpass)

Obs=y;
[nspecies, T]=size(Obs);
pstates=Posterior.gamma;
sequence=Sequence.sequence;
K=size(pstates,1);
time_bins=1:Lowpass:Lowpass*T;
ind=1:T;
if any(strcmp(fieldnames(Options),'lowpass'))
    lowpass=Options.lowpass;
end
ind=ind(1:1:T);
colors=colormap(parula(K+1));


if any(strcmp(fieldnames(Options),'id'))
    species_sort=Options.id;
end

subplot(2,1,1)

% 1st.

if 0
    % VITERBI
    hiddenPlot=Posterior.viterbiPath;
    hiddenPlot=hiddenPlot(ind);
    for k=1:K
      ndx=find(hiddenPlot==k);
      text(ndx, 0.5*ones(1,numel(ndx)),num2str(k), 'color', colors(k,:));
    end
end
% POSTERIOR
for k=1:K
    % prob
    h(k)=plot(1:T,pstates(k,1:T));
    hold on;
    set(h(k),'linewidth',2,'color',colors(k,:));
end   
%v=119*ones(1,101);
v= 28*ones(1,101);
plot(v,[0:.01:1],'linewidth',3,'color','r');
%v=125*ones(1,101);
v=42*ones(1,101);
plot(v,[0:.01:1],'linewidth',3,'color','r');
legend(h,strread(num2str(1:K),'%s'));

if 0
%     % shades
%     a=zeros(1,numel(ind));
%     for k=1:K
%         b=zeros(1,numel(ind));
%         indpstates=(pstates(k,ind)>Detection);
%         b(indpstates)=1;
%         [~,~]=jbfill(time_bins(ind),...
%             b,a,colors(k),0,0,0.2);
%     end
else
  if ~isempty(sequence)
        a=zeros(1,numel(ind));
        for k=1:size(sequence,2)
            b=zeros(1,numel(ind));
            indpstates=sequence(1,k):sequence(2,k);%(pstates(k,ind)>Detection);
            b(indpstates)=1;
            [~,~]=jbfill(time_bins(ind),...
                b,a,colors(sequence(4,k),:),0,0,0.1);
            hold on;
        end
        
    end
end
xlim([time_bins(ind(1)) time_bins(ind(end))]);
% 2nd.
%subplot(2,1,2);
% nspecies

subplot(2,1,2)
MaxObs=max(max(Obs));
Obs(y==0)=MaxObs+100;
Obs=log(Obs);
%imagesc(Obs);
%colordata(end,:) = [0.7 0.7 0.7];
xlabel('Time')
ylabel('species');
colors = colormap(parula(nspecies));
 for n=1:nspecies
     temp=Obs(n,ind);
     MAX=max(temp);
     MIN=min(temp);
     temp=(temp-MIN)/(MAX-MIN);
     hh=plot(time_bins(ind),temp);
     if(species_sort(n)<4)
        set(hh,'linewidth',1,'linestyle','-','color',colors(nspecies-n+1,:),'linewidth',3);
     else
         set(hh,'linewidth',1,'linestyle','-','color',colors(nspecies-n+1,:),'linewidth',1);
     end
     hold on
     v= 28*ones(1,101);
    plot(v,[0:.01:1],'linewidth',3,'color','r');

    v=42*ones(1,101);
    plot(v,[0:.01:1],'linewidth',3,'color','r');
             
  end
xlim([time_bins(ind(1)) time_bins(ind(end))]);
xlabel('time')
ylabel('species counts')
