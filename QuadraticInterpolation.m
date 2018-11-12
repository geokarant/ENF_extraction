function [ delta ] = QuadraticInterpolation( power_vector, index,fS )





for uu=1:numel(index);

alpha(1,uu) = ((power_vector((index(uu) - 1),uu)));
beta(1,uu) = ((power_vector((index(uu)),uu)));
gamma(1,uu) =((power_vector((index(uu) + 1),uu)));
delta1(1,uu) = (0.5*(alpha(1,uu) - gamma(1,uu)))./(alpha(1,uu) - 2*beta(1,uu) + gamma(1,uu));

delta(1,uu)=delta1(1,uu).*(fS(index(1,uu)+1)-fS(index(1,uu)))';


end
