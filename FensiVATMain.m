clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%GM1 data%%%%%%%%%%%%%%%%%%%%%%%
% dimension=1000;
%
% [data_matrix] = CS_data_generate_Punit(-6,1,25000,dimension);
% data_matrix_with_lables_1=[data_matrix zeros(25000,1)+1];
%
% [data_matrix] = CS_data_generate_Punit(0,2,50000,dimension);
% data_matrix=[data_matrix zeros(50000,1)+2];
% data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];
%
% [data_matrix] = CS_data_generate_Punit(6,3,25000,dimension);
% data_matrix=[data_matrix zeros(25000,1)+3];
% data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];
%
% Data=data_matrix_with_lables_1(:,1:end-1);
% Labels=data_matrix_with_lables_1(:,end);
%
% cp=10; ns=200;
% Reduced_Dim=50;
% Reduced_Dim_E=30;
% %%%%%%%%%GM2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimension=1000;
[data_matrix] = CS_data_generate_Punit(-2,1,25000,dimension);
data_matrix_with_lables_1=[data_matrix zeros(25000,1)+1];

[data_matrix] = CS_data_generate_Punit(0,2,50000,dimension);
data_matrix=[data_matrix zeros(50000,1)+2];
data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];

[data_matrix] = CS_data_generate_Punit(2,3,25000,dimension);
data_matrix=[data_matrix zeros(25000,1)+3];
data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];

Data=data_matrix_with_lables_1(:,1:end-1);
Labels=data_matrix_with_lables_1(:,end);


cp=20; ns=200;

Reduced_Dim=50; %%usually it needs to be estimated from JL lemma bound
% Reduced_Dim_E= Reduced_Dim; %%usually it needs to be estimated from JL lemma bound


n=size(Data,1);
INDX=randperm(n);
samples=n;
tt=INDX(1:samples);
Upspace_Mat=Data(tt,:);
lables= Labels(tt,:);
NoofK=length(unique(lables)); %%number of K from data.
clear tt INDX Data Labels %NormalizedData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DWNCPY=5; %%Number of RPs
N=1; %%Number of trials
Pi=zeros(n,1); %predicted labels
for i=1:N
    tic
    
    %%random projection
    T=rand(size(Upspace_Mat,2),Reduced_Dim);
    
    % %  First type of RP
    %     T(T<0.5)=-1;
    %     T(T>=0.5)=1;
    %      Downspace_Mat=Upspace_Mat*T/sqrt(Reduced_Dim);
    
    %   Second type of RP
    T(T<(1/3))=-sqrt(3);
    T(T>=(2/3))= sqrt(3);
    T(T<(2/3)&T>=(1/3))=sqrt(3);
    Downspace_Mat=Upspace_Mat*T/sqrt(Reduced_Dim);
    
    %third type or RP
    % Downspace_Mat = cos(.1/n*Upspace_Mat*(T));
    
    %% MMRS sampling
    smp= MMRS(Downspace_Mat, cp, ns );
    Samplingtime=toc;
    clear T;
    SampleData_UpSpace=Upspace_Mat(smp,:);
    SampleData_DWNSpace=Downspace_Mat(smp,:);
    Ub=ones(length(smp),length(smp));
    P=0; %%to be aggregated distance matrix
    tic
    for j=1:DWNCPY
        Reduced_Dim_EN= Reduced_Dim;
        %         %%multple RPs to get multiple distance matrix and cumulative aggregate them
        
        %% Alternative way to generate RP
        %         T=rand(size(SampleData_UpSpace,2),Reduced_Dim_E);
        %         T(T<0.5)=-0.5;
        %         T(T>=0.5)=0.5;
        %         downspace_Mat=SampleData_UpSpace*T/sqrt(Reduced_Dim_E);
        T=rand(size(Upspace_Mat,2),Reduced_Dim_EN);
        T(T<(1/3))=-sqrt(3);
        T(T>=(2/3))= sqrt(3);
        T(T<(2/3)&T>=(1/3))=sqrt(3);
        downspace_Mat=SampleData_UpSpace*T/sqrt(Reduced_Dim_EN);
        rs = distance2(downspace_Mat,downspace_Mat);
        Wrb= NormalizeRowsW(rs,1);
        Wrb=(Wrb+Wrb')/2;
        P=P+Wrb;
    end
    CumulativeENTime=toc;
    [labels_smp,Pi,ClusteringTime] = ClusteringMethods(Reduced_Dim,n,smp,Downspace_Mat,SampleData_DWNSpace,P,NoofK,'iVAT_Ensemble',Upspace_Mat,'DWNSPCTEST'); %%DWNSPC, DWNSPC_EN and DWNSPCTEST
    TotalTime_FastiVAT_EN=Samplingtime+CumulativeENTime+ClusteringTime;
    PA_New_Full_FastiVAT_EN= ClusterRelabellingPA(samples,NoofK,Pi,lables)
    
end
% Average_PA_New_Full_FastiVAT_EN=mean(PA_New_Full_FastiVAT_EN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
