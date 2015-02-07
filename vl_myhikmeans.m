function mytree = vl_myhikmeans(K, L,EnsembleWeight)
initialization
%mytree0=vl_mymytree( K, L); % saved in mytree0 to save memory
%mytree = vl_myfillhikmeans(Ensemble,K, L);
mytree = vl_myfillhikmeans_C(Ensemble,EnsembleWeight,K, L);
mytree.K=K;
mytree.depth=L;

if L==4
mytree.centers=int32(mytree.centers);
for i=1:K,
    mytree.sub(i).centers=int32(mytree.sub(i).centers);
    for j=1:K,
        mytree.sub(i).sub(j).centers=int32(mytree.sub(i).sub(j).centers);
        for k=1:K,
            mytree.sub(i).sub(j).sub(k).centers=int32(mytree.sub(i).sub(j).sub(k).centers);
        end
    end
end
end
if L==6
    mytree.centers=int32(mytree.centers);
    for i=1:K,
        mytree.sub(i).centers=int32(mytree.sub(i).centers);
        for j=1:K,
            mytree.sub(i).sub(j).centers=int32(mytree.sub(i).sub(j).centers);
            for k=1:K,
                mytree.sub(i).sub(j).sub(k).centers=int32(mytree.sub(i).sub(j).sub(k).centers);
                for k1=1:K,
                    mytree.sub(i).sub(j).sub(k).sub(k1).centers=int32(mytree.sub(i).sub(j).sub(k).sub(k1).centers);
                    for k2=1:K,
                        mytree.sub(i).sub(j).sub(k).sub(k1).sub(k2).centers=int32(mytree.sub(i).sub(j).sub(k).sub(k1).sub(k2).centers);
                    end
                end
            end
        end
    end
end

%save mytree mytree
