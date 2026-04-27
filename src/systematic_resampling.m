function ind = systematic_resampling(w)

N = length(w);

u0 = rand / N;
u = u0 + (0:N-1)'/N;

cw = cumsum(w);
ind = zeros(N,1);

i = 1;
j = 1;

while i <= N
    if u(i) <= cw(j)
        ind(i) = j;
        i = i + 1;
    else
        j = j + 1;
    end
end

end