function [chain, acc] = rto_metropolize(logrhos)

logrhos = logrhos(:)';
n   = length(logrhos);
ind = randperm(n);

tmp = logrhos(ind);
acc = (tmp(2:end) - tmp(1:end-1)) < log(rand(1,n-1));

chain = ones(1,n);
for i = 2:n
    if acc(i-1) 
        chain(i) = ind(i);
    else
        chain(i) = chain(i-1);
    end
end

end