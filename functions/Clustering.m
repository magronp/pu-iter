
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genom]=Clustering(B,G,M0,dBmax,beta,Fs,UseHierarchical)

if(nargin<7) UseHierarchical=false; end
if(nargin<6) Fs=44100; end
if(nargin<5) beta=1; end
if(nargin<4) dBmax=30; end
if(nargin<3) M=2; end

if(~UseHierarchical)
    [genom]=SingleStepClustering(B,M0,dBmax,beta,Fs);
else
    %%% algorithm for hierarchical clustering
    genom=ones(1,size(B,2));
    for M=2:M0
        % eval energies of separated clusters
        E=zeros(1,M-1);
        for m=1:length(E)
            mask=(genom==m);
            E(m)=sum(sum((B(:,mask)*G(mask,:)).^2));
        end

        % find a cluster with high energy and more than 1 element
        mask = 1;
        val  = 1;
        while(and(sum(mask)<2,val>0))
            [val,idx] = max(E);
            E(idx)    = 0;
            % separate the cluster with maximum energy further -> see also DAFx09
            mask      = (genom==idx);
        end
        if(sum(mask)<2)
            % error: no cluster with high energy and more than one element
            % found
            % -> no further subclustering
            subgenom=1;
        elseif(sum(mask)==2)
            % trivial partition of two signals in two clusters
            subgenom=1:2;
        else
            % further division in two subcluster
            subgenom=SingleStepClustering(B(:,mask),2,dBmax,beta,Fs);
        end
        % merge former clustering and new clustering
        step2=0;
        for step1=1:size(B,2)
            if(mask(step1))
                step2=step2+1;
                if(subgenom(step2)==2)
                    genom(step1)=M;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [genom]=SingleStepClustering(B,M,dBmax,beta,Fs)
R         = CreateFilterBank(2*size(B,1)-2,20,Fs);
tmp1_B    = R*(B.^2);
tmp2_B    = tmp1_B/max(tmp1_B(:))*(10^(dBmax/20)-1);
tmp3_B    = 20*log10(tmp2_B+1);

[Bcluster,Gcluster] = NMFcluster(tmp3_B,M,beta);
[val,genom]=max(Gcluster);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluates mel-filterbank
function [MelFiltMat]=CreateFilterBank(fftlen,N,Fs)
M      = round(fftlen/2)+1;
melmax = f2m(round(Fs/2));
s      = m2f((0:N-1)*melmax/(N-1));
k      = fftlen*s/Fs+1;
d      = diff(k);
c      = floor(k);
F      = zeros(N,M);
F(1,:) = (k(2)-(1:M))/d(1);
for iter=2:N-1
    F(iter,:)=[((1:c(iter))-k(iter-1))/d(iter-1),(k(iter+1)-(c(iter)+1:M))/d(iter)];
end
F(N,:)     = ((1:M)-k(N-1))/d(N-1);
MelFiltMat = max(F,0);

%%% avoids frequency spreading for low frequencies
for col=1:size(MelFiltMat,1)
    MelFiltMat(col,col)=MelFiltMat(col,col)+sum(MelFiltMat(col+1:end,col));
    MelFiltMat(col+1:end,col)=0;
end

function [B,G]=NMFcluster(X,I,beta)
%%% NMF with beta divergence
% random initialization
B = abs(randn(size(X,1),I));
G = abs(randn(I,size(X,2)));
for iter=1:100
    Xest=B*G;G=G.*(B'*((Xest.^(beta-2)).*X)+eps)./(B'*(Xest.^(beta-1))+eps); % update G
    Xest=B*G;B=B.*(((Xest.^(beta-2)).*X)*G'+eps)./((Xest.^(beta-1))*G'+eps); % update G
    
    % normalize energy between B and G
    for i=1:size(B,2)
        a=sqrt(sqrt(sum(G(i,:).^2))/sqrt(sum(B(:,i).^2)));
        B(:,i) = B(:,i)*a;
        G(i,:) = G(i,:)/a;        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=f2m(x)
y=log(1+x./700);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=m2f(x)
y=700*exp(x)-700;