   function [PA,cluster_matrix_mod] = ClusterRelabellingPA(samples,NoofK,PredictedLabel,lables)
    samples=length(lables);
    cluster_matrix_mod=zeros(1,samples);
    length_partition=zeros(1,NoofK);
    for i2=1:NoofK
        length_partition(i2)=length(find(PredictedLabel==i2));
    end
    [length_partition_sort,length_partition_sort_idx]=sort(length_partition,'descend');
    index_remaining=1:NoofK;
    
    for i2=1:NoofK
        original_idx=length_partition_sort_idx(i2);
        partition=find(PredictedLabel==original_idx);
        proposed_idx=mode(lables(partition));
        if(sum(index_remaining==proposed_idx)~=0)
            proposed_idx;
            cluster_matrix_mod(PredictedLabel==original_idx)=proposed_idx;
        else
            index_remaining(1);
            cluster_matrix_mod(PredictedLabel==original_idx)=index_remaining(1);
        end
        index_remaining(index_remaining==proposed_idx)=[];
    end
    
    PA=(samples-length(find((lables-(cluster_matrix_mod)'~=0))))/samples*100;