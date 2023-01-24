function [Kij2,e] = edge_selection(theta,SC,m,type)

if ~exist('type','var')
    disp('No edge selection paradigm was given: most synchronized will be used by default')
    type = 'Most';
end

n = size(SC,1);
N = size(SC,1);
inds = logical(triu(ones(N),1));
edges_flat = SC(inds);
edges_exist = find(edges_flat>0);

if type == 'Most'
    theta_i = theta * ones(1,n);                        % (nxn) matrix with same theta values in each row
    theta_j = ones(n,1) * theta';                       % (nxn) matrix with same theta values in each col
    r_i = exp(1i*theta_i);
    r_j = exp(1i*theta_j);
    r_ij = abs((r_i + r_j))./2; %take order parameter of all node pairs
    q = (SC>0).*r_ij;
    [~,I] = maxk(q,m);
    Kij2 = zeros(n);
    for ndx = 1:n
        maxkind = I(:,ndx);
        Kij2(maxkind,ndx) = SC(maxkind,ndx);
    end
    Kij2 = max(Kij2,Kij2');
    e_tmp = Kij2(inds)>0;
    e = e_tmp(edges_exist);

elseif type == 'Random'
    q = SC>0;
    z = 1:n;
    Kij2 = zeros(n);
    for ndx = 1:n
        I = q(:,ndx);
        if sum(I)<m
            inds = I;
        else
            inds = datasample(z(I),m,'Replace',false);
        end
        Kij2(inds,ndx) = SC(inds,ndx);
    end
    Kij2 = max(Kij2,Kij2');   
    e_tmp = Kij2(inds)>0;
    e = e_tmp(edges_exist);
    
elseif type == 'Least'
    theta_i = theta * ones(1,n);                        % (nxn) matrix with same theta values in each row
    theta_j = ones(n,1) * theta';                       % (nxn) matrix with same theta values in each col
    r_i = exp(1i*theta_i);
    r_j = exp(1i*theta_j);
    r_ij = abs((r_i + r_j))./2; %take order parameter of all node pairs
    q = (SC>0).*r_ij;
    q(SC==0) = inf;
    [~,I] = mink(q,m);
    Kij2 = zeros(n);
    for ndx = 1:n
        minkind = I(:,ndx);
        Kij2(minkind,ndx) = SC(minkind,ndx);
    end
    Kij2 = max(Kij2,Kij2');
    e_tmp = Kij2(inds)>0;
    e = e_tmp(edges_exist);  

else 
    disp('Type chosen is not an implemented one.')
end
    
    