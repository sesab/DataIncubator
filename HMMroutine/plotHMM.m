function plotHMM(y,Posterior,Sequence,Options,Lowpass)

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

% 1st.
subplot(5,1,1);
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
v=119*ones(1,101); %PATIENT B
%v=65*ones(1,101); %PATIENT A
plot(v,[0:.01:1],'linewidth',3,'color','r');
v=125*ones(1,101);%PATIENT B
%v=120*ones(1,101);%PATIENT A
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
                b,a,colors(sequence(4,k),:),0,0,0.2);
            hold on;
        end
        
    end
end
xlim([time_bins(ind(1)) time_bins(ind(end))]);

% 2nd.
subplot(5,1,[2 4]);
% nspecies

colormap default;
MaxObs=max(max(Obs));
Obs(y==0)=MaxObs+1000;
Obs=log(Obs);
imagesc(Obs);
hold on
l=species_sort(1)*ones(1,T);
plot([1:T],l(1,:),'k')
colormap default;
colordata = colormap ;
colordata(end,:) = [0.5 0.5 0.5];
colormap(colordata);
% colorbar;
xlabel('Time')
ylabel('species');
% for n=1:nspecies
%     temp=Obs(n,ind);
%     MAX=max(temp);
%     MIN=min(temp);
%     temp=(temp-MIN)/(MAX-MIN);
%     hh=plot(time_bins(ind),temp);
%     set(hh,'linewidth',1,'linestyle','--','color',colormap(n,:));
% end
xlim([time_bins(ind(1)) time_bins(ind(end))]);
hold off;






