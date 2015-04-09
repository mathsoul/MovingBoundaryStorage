    Policy = PolicyGenerator(result{3}{11},result{4}{11},NumX,NumQ);
    for i = 1 : NumQ
        for j = 1:NumX
            if(Policy(i,j)== 1)
                plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'o')
            elseif(Policy(i,j) == 2)
                plot(Qmin + (i-1)*(Qmax-Qmin)/(NumQ-1),Xmin + (j-1)*(Xmax-Xmin)/(NumX-1),'*')
            end        
        end
    end
%     
    plot([Qmin,Qmax],[alpha,alpha])

    axis([Qmin,Qmax,Xmin,Xmax])