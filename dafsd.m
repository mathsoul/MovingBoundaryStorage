for i = 1:Counter -1
    if sum(sum(Valuefunction{i+1}>Valuefunction{i}-0.1)) < 21*21
        disp(i)
    else
        mean(mean(Valuefunction{i+1}-Valuefunction{i}))
    end
end