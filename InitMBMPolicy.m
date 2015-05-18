function initial_policy = InitMBMPolicy()
    global NumQ NumX NumS
    initial_policy = zeros(NumQ,NumX,NumS);
    initial_policy(2:NumQ,NumX,:) = 2;
end