function [ rec_pos_trans2 ] = rotate_match_func( rec_pos_trans, true_pos_trans )
%This function finds the rotation (and potential reflection) to best match
%the reconstructed position to the true position. 

%The input are centered reconstructed positions (rec_pos_trans) and
%centered true positions (true_pos_trans). The output is the rotated
%reconstructed positions (rec_pos_trans2).

%First do rough grid search of possible rotations and reflections
dthetas=linspace(0,7/4*pi,7);
dsis=linspace(0,7/4*pi,7);
dphis=linspace(0,7/4*pi,7);
sse_min=Inf;
fmin=Inf;
for i=1:8
    %Test all possible reflection combinations
    temp_rec_pos_trans=rec_pos_trans;
    if i<5
        temp_rec_pos_trans(:,1)=-rec_pos_trans(:,1); %reflect across axis 1
    end
    if mod(i,2)==0
        temp_rec_pos_trans(:,2)=-rec_pos_trans(:,2); %reflect across axis 2
    end
    if i==3 || i==4 || i==7 || i==8
        temp_rec_pos_trans(:,3)=-rec_pos_trans(:,3); %reflect across axis 3
    end
    %Test all rotations given the current reflection
    for j=1:7
        for k=1:7
            for l=1:7
                dtheta=dthetas(j);
                dsi=dsis(k);
                dphi=dphis(l);
                sse=rotate_func([dtheta,dsi,dphi],temp_rec_pos_trans,true_pos_trans); %Get the sum-squared error or the rotation
                if sse<sse_min %If this is the best rotation for this reflection, use these parameters as initial parameters for a more thorough search
                    sse_min=sse;
                    init=[dtheta, dsi, dphi];
                end
            end
        end
    end
    %Find the best rotation (using fmincon instead of a gridsearch) for the
    %given reflection (using the gridsearch optimal parameters as the initial
    %parameters)
    [q,fval]=fmincon(@(x) rotate_func(x,temp_rec_pos_trans,true_pos_trans),init,[],[],[],[],[0 0 0],[2*pi 2*pi 2*pi]);
    %Now test whether this best fit (for this given reflection) is better
    %than the fits with other reflections. Save the parameters if this one
    %is the best
    if fval<fmin
        fmin=fval; %best sum squared error
        minVal=i; %best reflection param
        minQ=q; %best rotation param
    end
end

%Based on the optimal values found above, do the reflection and rotation

%Do the optimal reflection
if minVal<5
    rec_pos_trans(:,1)=-rec_pos_trans(:,1);
end
if mod(minVal,2)==0
    rec_pos_trans(:,2)=-rec_pos_trans(:,2);
end
if minVal==3 || minVal==4 || minVal==7 || minVal==8
    rec_pos_trans(:,3)=-rec_pos_trans(:,3);
end

%Do the optimal rotation
q=minQ;
theta=q(1);
psi=q(2);
phi=q(3);

%This is the rotation matrix in cartesian coordinates
R=[cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi);...
    -cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);...
    sin(theta) -sin(phi)*cos(theta) cos(phi)*cos(theta)];

rec_pos_trans2=rec_pos_trans*R;


end

