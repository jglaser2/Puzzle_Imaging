function [ rec_pos ] = sdm_func( W0 )
%This function does the SDM (sparse diffusion maps) algorithm to reduce dimensions

%Input similarity matrix (W0), and output 3-dimensional reconstructed
%positions (rec_pos)
%See manuscript (Table 1) for further description of all steps
%Steps listed below correspond to the manuscript

%Step 2
W0_norm=W0./repmat(sum(W0)',1,size(W0,1));
W=W0_norm;
%Step 3
opts.tol = 1e-7;
[v,e]=eigs(W,4,'lm',opts); %v are eigenvectors; e are eigenvales
rec_pos=v(:,2:4);
%Step 4
de=diag(e); %extract eigenvalues
rec_pos(:,1)=rec_pos(:,1)*de(2);
rec_pos(:,2)=rec_pos(:,2)*de(3);
rec_pos(:,3)=rec_pos(:,3)*de(4);

end

