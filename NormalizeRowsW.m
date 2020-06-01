function Wno= NormalizeRowsW(W,method)
switch method
    case 1
        for k=1:size(W,1)
            if sum(W(k,:))~0
                Wno(k,:)=W(k,:)./sum(W(k,:));
            else
                Wno(k,:)=  W(k,:);
            end
        end
    case 2
        
        for k=1:size(W,2)
            if sum(W(:,k))~0
                Wno(k,:)=W(:,k)./sum(W(:,k));
            else
                Wno(k,:)=  W(:,k);
            end
        end
end
end