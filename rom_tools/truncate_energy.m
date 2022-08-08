

function ind = truncate_energy( d, tol )

cd  = cumsum(d)/sum(d);
ind = sum( cd < 1-tol );

end