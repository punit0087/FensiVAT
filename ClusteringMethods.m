function [labels_smp,Pi,ClusteringTime] = ClusteringMethods(Reduced_Dim,n,smp,Downspace_Mat,DD,SampleData,NoofK,ClusteringMethod,Upspace_Mat,NPRMethod)
Pi=zeros(n,1);
kNN=3;% 3-5 works good
Pi_kNN=zeros(n,kNN);

switch ClusteringMethod
    case 'iVAT' %%corresponds to clusiVAT
        tic
        rs = distance2(SampleData,SampleData);
        [rv,C,I,ri,cut]=VAT(rs);
        [RiV,RV,reordering_mat]=iVAT(rv,1);
        
        %%% See Number of K from Image and Use Cut Threshold accordingly
        [cuts,ind]=sort(cut,'descend');
        ind=sort(ind(1:NoofK-1));
        
        Pi(smp(I(1:ind(1)-1)))=1;
        Pi(smp(I(ind(end):end)))=NoofK;
        for k=2:NoofK-1,
            Pi(smp(I(ind(k-1):ind(k)-1)))=k;
        end;
        SampleClusTime=toc;
               figure; imagesc(RiV); colormap(gray); axis image; axis off;
        labels_smp=Pi(smp);
        
    case 'iVAT_Ensemble' %%correspond to FensiVAT
        P=SampleData;
        tic
        [rv,C,I,ri,cut]=VAT(P);
        [RiV,RV,reordering_mat]=iVAT(rv,1);
        
        %%% See Number of K from Image and Use Cut Threshold accordingly
        [cuts,ind]=sort(cut,'descend');
        ind=sort(ind(1:NoofK-1));
        
        Pi(smp(I(1:ind(1)-1)))=1;
        Pi(smp(I(ind(end):end)))=NoofK;
        for k=2:NoofK-1
            Pi(smp(I(ind(k-1):ind(k)-1)))=k;
        end;
        SampleClusTime=toc;
        figure; imagesc(RiV); colormap(gray); axis image; axis off;
        labels_smp=Pi(smp);      
end

switch NPRMethod %%either to NPR in single downspace or voting from multiple RPs
    case 'DWNSPC' %%single downspace
        tic
        nsmp=setdiff(1:n,smp);
        r=distance2(Downspace_Mat(nsmp,:),DD);
        [~,s]=min(r,[],2);
        Pi(nsmp)=Pi(smp(s));
        FullDatalabelsTime=toc;
        
    case 'DWNSPC_EN' %%based on k-nearest NPRS (labels)
        tic
        nsmp=setdiff(1:n,smp);
        r=distance2(Downspace_Mat(nsmp,:),DD);
        [~,temp2]=sort(r,2);
        temp3=Pi(smp(temp2));
        Pi(nsmp)=mode(temp3(:,1:kNN),2);
        FullDatalabelsTime=toc;
    case 'UPSPC' %%NPR for clusiVAT
        tic
        nsmp=setdiff(1:n,smp);
        r=distance2(Upspace_Mat(nsmp,:),Upspace_Mat(smp,:));
        [~,s]=min(r,[],2);
        Pi(nsmp)=Pi(smp(s));
        FullDatalabelsTime=toc;
    case 'DWNSPCTEST' %% voting on multiple (k) RPs 
        tic
        nsmp=setdiff(1:n,smp);
        for i=1:kNN
            T=rand(size(Upspace_Mat,2),Reduced_Dim);
            T(T<0.5)=-1;
            T(T>=0.5)=1;
            Downspace_Matt=Upspace_Mat*T/sqrt(Reduced_Dim);
            DDD=Downspace_Matt(smp,:);
            r=distance2(Downspace_Matt(nsmp,:),DDD);
            [~,s]=min(r,[],2);
            Pi(nsmp)=Pi(smp(s));
            Pi_kNN(:,i)=Pi;
        end
        Pi=mode(Pi_kNN,2);
        FullDatalabelsTime=toc;
end
ClusteringTime= SampleClusTime+FullDatalabelsTime;
end