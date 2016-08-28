function M=computeDistanceMatrixGrid(nn,mm);

M=zeros(nn*mm);

X=[kron(ones(1,mm),1:nn);kron(1:mm,ones(1,nn))];

M=sqrt(pairwiseDistance(X));

%M=M/sum(M(:));
% MMM=sort(M(:),'ascend');
% baseTemp=median(MMM(nn*mm+1:end));
% M=M/baseTemp;

% below is faster for very large matrices.
% % if nn*mm > 60^2,
% %     baseTemp=median(median(M(1:15:end,1:11:end)))
% %     M=M/baseTemp;
% % else
% %     baseTemp=median(M(:));
% %     M=M/baseTemp;
% % end