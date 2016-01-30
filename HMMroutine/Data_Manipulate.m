
fileload='../Data/Mouse/Honda_mouse_data.mat';
%fileload='Data/LDdata/patientB.txt';
Dataload=load(fileload);
seq=Dataload;

ntrials=7;
nspecies=13;


DATA=repmat(struct('data',[]),ntrials,nspecies);
  for n=1:nspecies
    for t=1:ntrials
        DATA(t,n).data=seq.input.counts{t}(:,n);
        %temp=seq.input.counts(t)(:,CumSpecies(n));
        %temp=zscore(temp);
        % DATA(t,n).data=temp;
    end
  end
  
  v=zeros(56,nspecies*5);

  
      for i=1:nspecies
          n1=1;
          for n=1:5
          
            if(n1==4)
                  n1=n1+1;
              end
          for t=1:56
             
              n1
              n
              k=5*(i-1)+n;
                
             v(t,k)=DATA(n1,i).data(t);
             
          end
           n1=n1+1;
      end
     
  end
  
  dlmwrite('Mouse.txt',v,' ');