function [xnew,low,upp,lam,mu] = EIIM(m,n,loop,xval,xmin,xmax,...
    xold1,xold2,df0dx,fval,dfdx,low,upp,p,move,mu)
asy1=0.7; asy2=0.5; asyincr=5/6; asydecr=0.7; albefa=0.1; een=ones(n,1);
if loop<2.5
    low=xval-asy1*(xmax-xmin); upp=xval+asy2*(xmax-xmin);
else
    factor=asydecr*een; factor((xval-xold1).*(xold1-xold2)>0)=1/asyincr;
    low=xval-factor.*(xold1-low); upp=xval+factor.*(upp-xold1);
end
alfa=max(low+albefa*(xval-low),xmin); beta=min(upp-albefa*(upp-xval),xmax);
x_mu=(beta+alfa)/2; sigma=(beta-alfa)/2;
u=(xval-x_mu)./sigma;
umax=min(1,u+2*move); umin=max(-1,u-2*move);
if m==0
    h=0; dhdx=0*dfdx; mu=0;
elseif m==1
    h=fval; dhdx=dfdx;
elseif m>=2
    rho_ks=50;
    S=sum(exp(rho_ks*fval));
    h=(1/rho_ks)*log(S);
    dhdx=1/S*sum(exp(rho_ks*fval).*dfdx,2);
end
df0du=df0dx.*sigma; dhdu=dhdx.*sigma;
if m~=0 && loop>1.5 && h>=0
    F=@(mu)(mu*(u'*dhdu)+u'*df0du+norm(df0du+mu*dhdu,p/(p-1)));
    [mu]=fsolve(F,mu,optimset('Display','off'));
end
if h<0, mu=0; end
lam=norm(df0du+mu*dhdu,p/(p-1));
unew=max(umin,min(umax,-norm(u,p)*sign(df0du+mu*dhdu).*...
    (abs(df0du+mu*dhdu)/lam).^(1/(p-1))));
xnew=sigma.*unew+x_mu;