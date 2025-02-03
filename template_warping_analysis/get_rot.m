function [R,t]=get_rot(A,B)
cA=mean(A); % CoM A
cB=mean(B); % CoM B
N=size(A,1);
H=(A-repmat(cA,N,1))'*(B-repmat(cB,N,1)); % centred matrix H
[U,S,V]=svd(H); % SVD of H
R = V*U'; % finds rotation matrix

if det(R)<0 % if det is -ve then V needs altering
    V(:,3)=V(:,3).*-1;
    R=V*U';
end
t=-R*cA'+cB'; % find the translation
end
