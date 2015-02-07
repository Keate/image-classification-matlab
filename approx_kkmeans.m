function label=approx_kkmeans(Krect,k,max_iter,indices)   
% Krect - m X N rectangular kernel matrix  
% k - number of clusters
% max_iter - maximum number of iterations
% indices - indices of the sampled data points

    [m,N]=size(Krect);
    Khat=Krect(1:m,indices);
    Krect=full(Krect);
    
    Khat_inv=pinv(full(Khat));
    T = Khat_inv*Krect;
    
    init_labels = single(ceil(rand(1,N)*k));
    label =init_labels;
    
    last=0;
    t=0;
    U=zeros(max(label),N,'single');
    while(any(label~=last) && t < max_iter)
        U=U-U;
        for k=1:length(label),
            U(label(k),k)=1;
        end
% 	    U=sparse(label,1:N,1,k,N,N);  % compute sparse membership
%         U=full(U);
	    U=full(bsxfun(@rdivide,U,sum(U,2)));
        
	    alpha=(T*U')';  % solve for cluster center weights
	    D=bsxfun(@plus,-2*alpha*Krect,diag(alpha*Khat*alpha')); % distance between cluster center and each object
	    last=label;
	    [temp,label]=min(D); % assign new labels
	        
	    t=t+1;
    end;
    if(t >= max_iter)
	    fprintf('Approx kk means - Max iter exceeded');
    end
end
