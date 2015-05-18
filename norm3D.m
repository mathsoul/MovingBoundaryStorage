function norm_3D = norm3D(policy_diff)
    global NumQ NumX NumS
    
    policy_diff = reshape(policy_diff,[NumQ*NumX*NumS,1]);
    norm_3D = norm(policy_diff);
end
