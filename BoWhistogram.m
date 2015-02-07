function hist_test=BoWhistogram(f,codewordnum,C,knorm)
% 
% 
%             index=zeros(1,size(f,1));
%             F=zeros(size(C));
            if 1==knorm,
%                 for i=1:size(f,1),
%                     F=repmat(f(i,:),codewordnum,1);
%                     D=abs(F-C);
%                     dis=sum(D,2);
%                     [mindis,index(i)]=min(dis);
%                 end
                [drop, index] = min(vl_alldist(C', double(f'),'l1'), [], 1) ;
            else
                [drop, index] = min(vl_alldist(C', double(f'),'l2'), [], 1) ;
            end
            
hist_test=hist(index,[1:codewordnum]);
hist_test=hist_test/(sum(hist_test)+eps);
