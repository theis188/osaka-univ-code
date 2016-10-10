function [Us,conc,bif,TimeTaken] = SimpleODESolverMatlab(NoSteps,InitialState,U0,Uf, ...
    Kvec,DVDX,DVDU,S)
%SimpleODESolver: Fixed Step No Correction Solver
options = odeset('Events',@event_function);
TimeIn = clock;
[t, conc,TE,XE] = ode15s(@dxFunc,[0:1/NoSteps:1],InitialState,options,Kvec,U0,Uf,DVDX,DVDU,S,TimeIn);

Us = repmat(U0,1,length(t)) + repmat(t',length(U0),1).*repmat(Uf-U0,1,length(t));
bif = t(end) ~= 1;

TimeTaken = etime(clock,TimeIn);
end

function dx = dxFunc(t,x,Kvec,Uini,Ufinal,DVDX,DVDU,S,~)
    U = Uini + t.*(Ufinal-Uini);
    XJac = S*DVDX(x,Kvec,1,U);
    dx = -(XJac)\S*DVDU(x,Kvec,1,U)*(Ufinal-Uini);
end

function [value,isterminal,direction] = event_function(t,x,Kvec,Uini,Ufinal,DVDX,~,S,TimeIn)
U = Uini + t.*(Ufinal-Uini);
XJac = S*DVDX(x,Kvec,1,U);
value = max(real(eig(XJac)));  % when value = 0, an event is triggered]
if max(real(eig(XJac)))>0 || any(x<0) %|| etime(clock,TimeIn)>30,
    value = 0;
end
isterminal = 1; % terminate after the first event
direction = 0;  % get all the zeros 
end