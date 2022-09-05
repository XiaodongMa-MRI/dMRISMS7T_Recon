function x = fnlCg_l1spirit(x0,params)
%--------------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F* W' *x - y||^2 + lambda1*|x|_1 + lambda2*TV(W'*x) 
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% (c) Michael Lustig 2007
%--------------------------------------------------------------------------
x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha; ,    beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

% copmute g0  = grad(Phi(x))

g0 = wGradient(x,params);

dx = -g0;


% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);
	f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, 0, params);
	t = t0;
    
        f1  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, t, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx);

	%--------- uncomment for debug purposes ------------------------	
    if params.dsp
    	disp(sprintf('%d   , obj: %f', k, f1*100));
    end
    %---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
		break;
	end

end


end



%--------------------------------------------------------------------------
function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)
% precalculates transforms to make line search cheap

X = zeros([params.N, params.num_chan]);
X(params.idx_mis) = x;

FTXFMtx = params.GOP * X;


dX = zeros([params.N, params.num_chan]);
dX(params.idx_mis) = dx;

FTXFMtdx = params.GOP * dX;


if params.TVWeight
    DXFMtx = params.TV * X;
    DXFMtdx = params.TV * dX;
else
    DXFMtx = 0;
    DXFMtdx = 0;
end

end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function [res, obj] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, t, params)
%calculates the objective function

obj = FTXFMtx + t*FTXFMtdx - params.data;
obj = obj(:)'*obj(:);

if params.TVWeight
    w = DXFMtx + t*DXFMtdx + params.grad;
    
    if params.L == 1
        % L1 over coils
        TV = abs(w); 
    else
        if params.L == 12
            % joint L21 over coils
            TV = sum(abs(w).^2, 3).^0.5; 
        end
    end
    
else
    TV = 0;
end


TV = sum(TV(:)) * params.TVWeight;

res = obj + TV;

end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function grad = wGradient(x,params)

gradTV = 0;

gradObj = gOBJ(x,params);


if params.TVWeight
    gradTV = gTV(x,params);
end


grad = gradObj + params.TVWeight * gradTV;

end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

X = zeros([params.N, params.num_chan]);
X(params.idx_mis) = x;

gradObj = params.GOP' * ((params.GOP * X) - params.data);

gradObj = 2*gradObj(params.idx_mis);

end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function grad = gTV(x,params)
% computes gradient of TV operator

X = zeros([params.N, params.num_chan]);
X(params.idx_mis) = x;

Dx = (params.TV * X) + params.grad;

if params.L == 1
    G = sign(Dx); 
else
    if params.L == 12
        DxL2 = sum(abs(Dx).^2, 3).^0.5;
        G = Dx ./ (eps + repmat(DxL2, [1,1,params.num_chan,1]));
    end
end

grd = params.TV' * G;

grad = grd(params.idx_mis);

end
%--------------------------------------------------------------------------





