% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [P,residual,fitcurve,forecastcurve,timevect2,initialguess,fval,flagg,fitDiagnostics]=fit_model(data1,params0,fixI0,numstartpoints,DT,flagX,forecastingperiod,fitopts)
%FIT_MODEL  Fit one growth model to one calibration window.
%
%   [P,residual,fitcurve,forecastcurve,timevect2,initialguess,fval,flagg,fitDiagnostics] = ...
%       fit_model(data1,params0,fixI0,numstartpoints,DT,flagX,forecastingperiod)
%   [...] = fit_model(..., fitopts)
%
%   P = [r p a K I0 alpha d]. The first eight outputs are as before.
%
%   Optional 8th input FITOPTS, a struct; every field may be omitted:
%     fitopts.logscale      true (default): the optimizer and the random start
%                           points work in log10 of r, K and alpha (the
%                           coordinates that span orders of magnitude).
%     fitopts.useParallel   false (default); true runs the MultiStart local
%                           solves on a parallel pool (Parallel Computing
%                           Toolbox), falling back to serial with one warning.
%     fitopts.optimOptions  struct of fmincon option overrides, for example
%                           struct('StepTolerance',1e-8).
%     fitopts.LB, fitopts.UB  1-by-7 bounds that replace the bounds derived
%                           from data1, and make params0 the seed as given.
%                           The runners pass the bounds of the fit to the
%                           observed data to every bootstrap refit (see below).
%
%   Optional 9th output FITDIAGNOSTICS: bounds, convergence flags, the share of
%   converged starts that agree with the best one (startAgreement) and which
%   estimated coordinates ended on a bound (atLB, atUB).
%
%   REVISED November 2026, following fit_model in QuantDiffForecast:
%
%   Optimizer coordinates. fmincon no longer receives z = [r p a K I0 alpha d].
%   It receives only the coordinates that are estimated (UB > LB), each mapped
%   to the unit interval, u = (t - LBt)./(UBt - LBt), with t = z or log10(z)
%   (see fitopts.logscale), and the objective divided by its value at the seed.
%   Reason: the sqp step tolerance is relative to the size of the whole
%   decision vector. With K ~ 1e4-1e9 in that vector, any step shorter than
%   about K*StepTolerance passed the test, so a local run could stop with a
%   positive exit flag after moving r or p by far less than their uncertainty.
%   In unit coordinates the tolerances mean what they say. Bounds, start
%   points, P and fval are all returned in the original units.
%
%   Start points. numstartpoints points cover the box (Latin hypercube when
%   lhsdesign is available, uniform otherwise) in the optimizer's coordinates,
%   so K is explored on the log scale instead of almost always near its upper
%   bound; up to three jittered copies of the seed and the seed itself are added.
%
%   Failed solves. growthObjective returns Inf for a trial that cannot be
%   evaluated; localBound maps it above every valid value (was: 10^10).
%
%   One evaluation. fval, fitcurve and residual come from the same call of
%   growthObjective at the returned P.
%
%   Bounds of bootstrap refits. The bounds of r and K are derived from the data
%   (r within a factor 20 of the slope of log incidence over the first 7 points,
%   K at or above the cumulative count). Recomputed for every bootstrap
%   replicate, the r interval is driven by the noise in 7 small counts, and for
%   some replicates it excluded the fitted r altogether (about 4% of replicates
%   of the shipped series curve-GLM...dist1-1-factor1-1, where r was then pinned
%   as low as 0.017 against a fit of 0.17). With fitopts.LB/UB all refits share
%   the box of the observed-data fit, as in QuantDiffForecast, and start from
%   the supplied params0 (the best fit) instead of a re-derived r.
%
%   Binding interval for r. For the fit to the observed data the same data-derived
%   interval can exclude the right growth rate (1 of 60 simulated GLM series in
%   the tests: r in [4e-5, 0.017] for a true r of 0.18). When the fit ends on that
%   bound it is repeated once with r in [1e-3, 10]; see the block after the final
%   evaluation and fitDiagnostics.rBoundsWidened.
%
%   Logistic model with fixI0 = 0: the unused exponent a had bounds [1 10] and
%   was therefore estimated; it is now fixed at 1 like in the fixI0 = 1 case.

global flag1 method1 timevect ydata

if nargin < 8
    fitopts = [];
end

flag1=flagX;

if isempty(method1)
    error('GrowthPredict:fit_model:NoMethod', ...
        'The global method1 is empty. It is set inside the options file (options_fit / options_forecast).');
end

if any(~isfinite(data1(:)))
    error('GrowthPredict:fit_model:InvalidData', ...
        'fit_model received non-finite data (NaN/Inf).');
end

timevect=data1(:,1)*DT;

r=params0(1);
p=params0(2);
a=params0(3);
K=params0(4);
I0=params0(5);
alpha=params0(6);
d=params0(7);

I0=data1(1,2); % initial condition

z(1)=r;
z(2)=p;
z(3)=a;
z(4)=K;
z(5)=I0;
z(6)=alpha;
z(7)=d;


switch method1
    case 0
        LBe=[0 0];
        UBe=[0 0];
    case 1
        LBe=[0 0];
        UBe=[0 0];
    case 3
        LBe=[10^-8 1];
        UBe=[10^3 1];
    case 4
        LBe=[10^-8 1];
        UBe=[10^5 1];
    case 5
        LBe=[10^-8 0.6]; %d>=1
        UBe=[10^5 10^3];
    otherwise
        error('GrowthPredict:fit_model:UnsupportedMethod', ...
            'method1 = %g is not supported; use 0, 1, 3, 4 or 5.', method1);

end

w  = min(7, size(data1,1));                    % early window

if w < 2                                       % degenerate: too few points
    r0  = params0(1);
    rlb = 1e-3;  rub = 10;
else
    y  = max(data1(1:w,2), 0.5);               % guard zeros
    b  = polyfit(data1(1:w,1)*DT, log(y), 1);

    if b(1) <= 0
        r0  = params0(1);      % no growth signal: keep supplied guess
        rlb = 1e-3;  rub = 10; % wide default
    else
        r0  = b(1);
        rlb = r0/20; rub = r0*20;
    end
end

Cobs = sum(abs(data1(:,2)));
Klb  = max(20, Cobs);            % can't undershoot observed cumulative
Kub  = min(1e3*Cobs, 1e9);       % capped
Kub  = max(Kub, Klb);            % keep Klb <= Kub for very small series

if fixI0==1

    switch flag1

         case -1  %EXP
            LB=[rlb  1 1 0 I0 LBe];
            UB=[rub  1 1 0 I0 UBe];

        case 0   %GGM
            LB=[rlb  0.01 1 0 I0 LBe];
            UB=[rub  1 1 0 I0 UBe];

        case 1 % GLM
            LB=[rlb  0.01 1 Klb I0 LBe];
            UB=[rub  1 1 Kub I0 UBe];

        case 2 %GRM
            LB=[rlb  0.01 0 Klb I0 LBe];
            UB=[rub  1 10 Kub I0 UBe];

        case 3 %Logistic
            LB=[rlb  1 1 Klb I0 LBe];
            UB=[rub  1 1 Kub I0 UBe];

        case 4 % Richards
            LB=[rlb  1 0 Klb I0 LBe];
            UB=[rub  1 10 Kub I0 UBe];

        case 5 % Gompertz
            LB=[max(1e-4, params0(1)/10)  1 max(1e-6, params0(3)/10) 1 I0 LBe];
            UB=[params0(1)*10             1 params0(3)*10            1 I0 UBe];

        otherwise
            error('GrowthPredict:fit_model:UnknownModel', ...
                'flag1 = %g is not a supported growth model (-1 to 5).', flag1);

    end

else

    y1 = data1(1,2);

    I0lb=max(1, 0.2*y1);
    I0ub=5*max(1,y1);
    I0ub=min(I0ub, Klb);         % keep I0 <= K for saturating models
    I0lb=min(I0lb, I0ub);        % preserve I0lb <= I0ub

    switch flag1

        case -1
            LB=[rlb 1 1 0 I0lb LBe];
            UB=[rub  1 1 0 I0ub UBe];

        case 0
            LB=[rlb 0.01 1 0 I0lb LBe];
            UB=[rub  1 1 0 I0ub UBe];

        case 1
            LB=[rlb  0.01 1 Klb I0lb LBe];
            UB=[rub 1 1 Kub I0ub UBe];

        case 2
            LB=[rlb  0.01 0 Klb I0lb LBe];
            UB=[rub  1 10 Kub I0ub UBe];

        case 3 % Logistic: the exponent a is not part of the model, keep it fixed at 1
            LB=[rlb  1 1 Klb I0lb LBe];
            UB=[rub  1 1 Kub I0ub UBe];

        case 4
            LB=[rlb  1 0 Klb I0lb LBe];
            UB=[rub  1 10 Kub I0ub UBe];

        case 5 % Gompertz
            LB=[max(1e-4, params0(1)/10)  1 max(1e-6, params0(3)/10) 1 I0lb LBe];
            UB=[params0(1)*10             1 params0(3)*10            1 I0ub UBe];

        otherwise
            error('GrowthPredict:fit_model:UnknownModel', ...
                'flag1 = %g is not a supported growth model (-1 to 5).', flag1);

    end

end

if flag1 ~= 5
    z(1) = min(max(r0, LB(1)), UB(1));   % informed r guess
end

% Caller-supplied box (bootstrap refits): use it as is and keep params0 as the seed.
boundsSupplied = isstruct(fitopts) && isfield(fitopts,'LB') && isfield(fitopts,'UB') && ...
    ~isempty(fitopts.LB) && ~isempty(fitopts.UB);
if boundsSupplied
    if numel(fitopts.LB)~=numel(z) || numel(fitopts.UB)~=numel(z)
        error('GrowthPredict:fit_model:SuppliedBounds', ...
            'fitopts.LB and fitopts.UB must have %d entries [r p a K I0 alpha d].', numel(z));
    end
    LB=fitopts.LB;
    UB=fitopts.UB;
    z=params0;
end

LB=LB(:)'; UB=UB(:)'; z=z(:)';
nz=numel(z);

if numel(LB)~=nz || numel(UB)~=nz
    error('GrowthPredict:fit_model:Length', ...
        'Length mismatch: numel(z)=%d, numel(LB)=%d, numel(UB)=%d.', nz, numel(LB), numel(UB));
end
if any(LB>UB) || ~all(isfinite(LB)) || ~all(isfinite(UB))
    error('GrowthPredict:fit_model:Bounds', ...
        'The parameter bounds must be finite with LB <= UB (check the data and the initial guesses).');
end

z0 = min(max(z, LB), UB);                % clamp everything into bounds
mid=(LB+UB)/2;
z0(isnan(z0))=mid(isnan(z0));

ydata=data1(:,2);

% <============================================================================>
% <================= Optimizer coordinates: free parameters, unit box =========>
% <============================================================================>

useLog      = localFlag(fitopts,'logscale',true);
useParallel = localFlag(fitopts,'useParallel',false);

free  = find(UB > LB);
nFree = numel(free);

tf = false(1,nz);                        % coordinates searched in log10
if useLog
    logCandidates = [1 4 6];             % r, K, alpha
    tf(logCandidates) = (UB(logCandidates) > LB(logCandidates)) & (LB(logCandidates) > 0);
end

LBt = localToT(LB,tf);  UBt = localToT(UB,tf);  spanT = UBt-LBt;
toQ   = @(Z) localToQ(Z,tf,free,LBt,spanT);          % rows of z -> rows of u
fromQ = @(U) localFromQ(U,tf,free,LBt,spanT,z0);     % rows of u -> rows of z

% Objective in original units; no global variables, so it can run on workers.
tvecFit=timevect; yFit=ydata; flagFit=flag1; methodFit=method1;
f  = @(Z) growthObjective(Z,tvecFit,yFit,flagFit,methodFit);
fq = @(U) f(fromQ(U));

f0 = f(z0);
if isfinite(f0) && abs(f0) >= 1
    fscale = abs(f0);
else
    fscale = 1;
end
fqs = @(U) localBound(fq(U)/fscale);     % what fmincon minimizes

options=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'StepTolerance',1e-10, ...
    'FunctionTolerance',1e-10, ...
    'OptimalityTolerance',1e-8, ...
    'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-6, ...
    'ScaleProblem','none', ...           % scaling is done explicitly above
    'MaxFunctionEvaluations',20000, ...
    'MaxIterations',20000, ...
    'Display','off');

if isstruct(fitopts) && isfield(fitopts,'optimOptions') && isstruct(fitopts.optimOptions)
    keys=fieldnames(fitopts.optimOptions);
    for j=1:numel(keys)
        options=optimoptions(options,keys{j},fitopts.optimOptions.(keys{j}));
    end
end

% <============================================================================>
% <============================ Start points ==================================>
% <============================================================================>

nStarts = max(0, floor(numstartpoints));
nJit    = min(3, nStarts);
lhsOK   = exist('lhsdesign','file') == 2;

makeStarts = @() localFromT(localMakeStarts(nStarts,nJit,LBt,UBt,localToT(z0,tf),free,lhsOK), tf);

starts = makeStarts();
initialguess = starts;

persistent parallelUnavailable
ms = MultiStart('Display','off');
if useParallel && ~isequal(parallelUnavailable,true)
    try
        gcp;                                           % start or reuse the pool
        ms = MultiStart('Display','off','UseParallel',true);
    catch ME
        parallelUnavailable = true;
        warning('GrowthPredict:fit_model:Parallel', ...
            ['Parallel pool unavailable (%s); running MultiStart serially for the rest of ' ...
            'this session (clear fit_model to retry).'], ME.message);
    end
end

maxAttempts=3; % Initial attempt plus at most two retries.
attempt=0;
flagg=-Inf;
fval=Inf;
P=z0;
bestP=[]; bestFval=Inf;
startAgreement=NaN;

if nFree==0

    % Nothing to estimate: evaluate the fixed vector once.
    fval=f(z0); flagg=1; bestP=z0; bestFval=fval; startAgreement=1;

else

    problem = createOptimProblem('fmincon','objective',fqs,'x0',toQ(z0), ...
        'lb',zeros(1,nFree),'ub',ones(1,nFree),'options',options);

    % Retry with fresh start points when MultiStart reports failure, or when
    % its best point is a failed solve (runs that start where every solve
    % fails sit on a flat plateau and report convergence).
    while (flagg<=0 || ~isfinite(fval)) && attempt<maxAttempts

        attempt=attempt+1;

        if attempt>1
            starts=makeStarts();
            initialguess=[initialguess;starts]; %#ok<AGROW>
        end

        sp = CustomStartPointSet(toQ(starts));

        [Pq,fvalq,flagg,~,allmins] = run(ms,problem,sp);

        if isempty(Pq)
            fval=Inf;
        else
            P=fromQ(Pq);
            fval=localUnbound(fvalq)*fscale;
        end

        if ~isempty(allmins)
            fvs=arrayfun(@localUnbound,[allmins.Fval])*fscale;
            [fv,idx]=min(fvs);
            if isfinite(fv) && fv<bestFval
                bestFval=fv;
                bestP=fromQ(allmins(idx).X);
            end
            okRun=([allmins.Exitflag]>0) & isfinite(fvs);
            if any(okRun)
                startAgreement=sum(okRun & fvs<=bestFval+1e-4*max(1,abs(bestFval)))/sum(okRun);
            end
        end

    end

end

if flagg>0 && ~isfinite(fval)
    flagg=0;   % a "converged" failed solve is not a solution
end

if flagg<=0
    if ~isempty(bestP)
        warning('GrowthPredict:fit_model:MultiStart', ...
            'MultiStart did not converge cleanly after %d attempts; using the best local minimum found (fval = %g).', ...
            attempt,bestFval);
        P=bestP; fval=bestFval;
    else
        error('GrowthPredict:fit_model:MaxAttemptsReached', ...
            ['Fit failed after %d attempts (last exit flag: %g): no start point gave a valid objective value. ' ...
            'Common causes: a first data point equal to zero with fixI0 = 1, or data the model cannot follow.'], ...
            attempt,flagg);
    end
end

% Rescale-and-confirm pass. The objective was divided by |f(z0)|. When the seed
% is far from the optimum, the scaled objective at the optimum is tiny and
% fmincon's relative tolerances were much looser than intended (and the other
% way round when the seed could not be evaluated). If the scale is off by more
% than a factor of 100, run one more local solve from P with the scale reset.
rescale=struct('run',false,'scaleBefore',fscale,'scaleAfter',fscale,'exitflag',NaN,'improved',false);
fscaleNew=max(abs(fval),1);
if nFree>0 && isfinite(fval) && (fscaleNew<1e-2*fscale || fscaleNew>1e2*fscale)
    rescale.run=true; rescale.scaleAfter=fscaleNew;
    try
        fqsNew=@(U) localBound(fq(U)/fscaleNew);
        [qc,fc2,flc]=fmincon(fqsNew,toQ(P),[],[],[],[],zeros(1,nFree),ones(1,nFree),[],options);
        qc=min(max(qc(:)',0),1);
        fc2=localUnbound(fc2)*fscaleNew;
        rescale.exitflag=flc;
        if isfinite(fc2) && fc2<fval
            P=fromQ(qc); fval=fc2; rescale.improved=true;
        end
    catch
        % keep the MultiStart point
    end
    fscale=fscaleNew;
end

% Final evaluation: fval, fitcurve and residual from ONE call at the returned P.
% This call throws if the curve cannot be evaluated.
optimizerFval=fval;
[fval,fitcurve]=growthObjective(P,timevect,ydata,flag1,method1);

residual=fitcurve-ydata;

tolB=1e-3;
atLB=false(1,nz); atUB=false(1,nz);
atLB(free)=P(free)<=LB(free)+tolB*abs(LB(free));
atUB(free)=P(free)>=UB(free)-tolB*abs(UB(free));

% Escape from a binding data-derived interval for r. The interval [r0/20, 20*r0]
% comes from the slope of log incidence over the first 7 points; with small early
% counts that slope can be far too small and the interval then excludes the
% growth rate the rest of the window supports (the fit ends on the r bound with
% p and K on their bounds as well). When the fit lands on that bound, refit once
% from the current solution with r allowed in [1e-3, 10] (the default interval
% used when there is no growth signal) and keep the result if it is better. The
% widened bounds are returned in fitDiagnostics, so bootstrap refits use them too.
if ~boundsSupplied && flag1~=5 && (atLB(1) || atUB(1)) && (LB(1)>1e-3 || UB(1)<10)

    wideopts=fitopts;
    if ~isstruct(wideopts)
        wideopts=struct();
    end
    wideopts.LB=LB;
    wideopts.UB=UB;
    wideopts.LB(1)=min(LB(1),1e-3);
    wideopts.UB(1)=max(UB(1),10);

    [P2,residual2,fitcurve2,forecastcurve2,timevect22,initialguess2,fval2,flagg2,fitDiagnostics2]= ...
        fit_model(data1,P,fixI0,numstartpoints,DT,flagX,forecastingperiod,wideopts);

    if fval2<fval
        warning('GrowthPredict:fit_model:RBoundsWidened', ...
            ['The fit ended on the data-derived bound of r (%g in [%g, %g]); refitted with r in [%g, %g]: ' ...
            'r = %g, objective %g -> %g.'],P(1),LB(1),UB(1),wideopts.LB(1),wideopts.UB(1),P2(1),fval,fval2);
        P=P2; residual=residual2; fitcurve=fitcurve2; forecastcurve=forecastcurve2; timevect2=timevect22;
        initialguess=[initialguess;initialguess2]; fval=fval2; flagg=flagg2;
        fitDiagnostics=fitDiagnostics2;
        fitDiagnostics.rBoundsWidened=true;
        return
    end

    % the wider interval did not help: restore the state of this call
    timevect=data1(:,1)*DT;
    ydata=data1(:,2);

end

if nargout>8
    fitDiagnostics=struct('LB',LB,'UB',UB,'freeIndex',free,'nFree',nFree, ...
        'logscaleMask',tf,'seed',z0,'objectiveAtSeed',f0,'objectiveScale',fscale, ...
        'optimizerObjective',optimizerFval,'finalObjective',fval, ...
        'exitflag',flagg,'attempts',attempt,'startAgreement',startAgreement, ...
        'numberOfStarts',size(starts,1),'rescalePass',rescale, ...
        'atLB',atLB,'atUB',atUB,'boundTolerance',tolB,'boundsSupplied',boundsSupplied, ...
        'rBoundsWidened',false);
end

if forecastingperiod<1

    forecastcurve=fitcurve;
    timevect2=timevect;

else

    timevect2=(data1(1,1):data1(end,1)+forecastingperiod)*DT;

    forecastcurve=growthIncidenceCurve(P,timevect2,flag1);

end

end % fit_model


% <============================================================================>
% <============================ Local functions ===============================>
% <============================================================================>

function val=localFlag(fitopts,name,default)
val=default;
if isstruct(fitopts) && isfield(fitopts,name) && ~isempty(fitopts.(name))
    val=logical(fitopts.(name));
end
end

function starts=localMakeStarts(nBox,nJit,LBt,UBt,t0,free,lhsOK)
% Start points in the optimizer's (possibly log10) coordinates, one per row.
%   nBox points cover the box of the free coordinates,
%   nJit points are jittered copies of the seed (sd = 5% of each range),
%   the seed itself is always the last row.
nz=numel(t0);
nF=numel(free);
span=UBt-LBt;

S1=repmat(t0,nBox,1);
if nBox>0 && nF>0
    if lhsOK
        X=lhsdesign(nBox,nF,'criterion','maximin','iterations',50);
    else
        X=rand(nBox,nF);
    end
    S1(:,free)=repmat(LBt(free),nBox,1)+X.*repmat(span(free),nBox,1);
end

S2=repmat(t0,nJit,1);
if nJit>0 && nF>0
    S2(:,free)=S2(:,free)+0.05*randn(nJit,nF).*repmat(span(free),nJit,1);
    S2=min(max(S2,repmat(LBt,nJit,1)),repmat(UBt,nJit,1));
end

starts=[S1;S2;t0];
starts=unique(starts,'rows','stable');
starts=reshape(starts,[],nz);
end

function T=localToT(Z,tf)
T=Z;
if any(tf)
    T(:,tf)=log10(max(Z(:,tf),realmin));
end
end

function Z=localFromT(T,tf)
Z=T;
if any(tf)
    Z(:,tf)=10.^T(:,tf);
end
end

function U=localToQ(Z,tf,free,LBt,spanT)
% Rows of Z (original units) -> rows of U (unit box, free coordinates only).
T=localToT(Z,tf);
n=size(T,1);
U=(T(:,free)-repmat(LBt(free),n,1))./repmat(spanT(free),n,1);
U=min(max(U,0),1);
end

function Z=localFromQ(U,tf,free,LBt,spanT,z0)
% Rows of U (unit box, free coordinates) -> rows of Z (original units, full
% vector). Fixed coordinates are taken from the seed z0.
nF=numel(free);
if size(U,2)~=nF && size(U,1)==nF
    U=U.';                               % a single point passed as a column
end
n=size(U,1);
U=min(max(U,0),1);
T=repmat(localToT(z0,tf),n,1);
T(:,free)=repmat(LBt(free),n,1)+U.*repmat(spanT(free),n,1);
Z=localFromT(T,tf);
end

function g=localBound(v)
% Monotone, bounded transform of the scaled objective handed to fmincon.
%   v <= T       : g = v
%   v >  T       : g = T + log1p(v - T)   (continuous, <= T + 710)
%   v not finite : g = T + FAILGAP        (a failed solve)
% A failed trial is therefore worse than every finite objective value.
T=1e6; FAILGAP=1e3;
if ~isfinite(v) || ~isreal(v)
    g=T+FAILGAP;
elseif v>T
    g=T+log1p(v-T);
else
    g=v;
end
end

function v=localUnbound(g)
% Inverse of localBound; the failure value maps back to Inf.
T=1e6; FAILGAP=1e3;
if ~isfinite(g) || g>=T+FAILGAP
    v=Inf;
elseif g>T
    v=T+expm1(g-T);
else
    v=g;
end
end
