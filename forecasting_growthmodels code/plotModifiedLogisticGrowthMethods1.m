function objfunction=plotModifiedLogisticGrowthMethods1(z)
%PLOTMODIFIEDLOGISTICGROWTHMETHODS1  Legacy interface to the fitting objective.
%
%   Kept so that existing scripts calling this function keep working. It reads
%   the global variables set by fit_model and forwards to growthObjective,
%   which fit_model now calls directly. A trial that cannot be evaluated
%   returns Inf (formerly 10^10).

global flag1 method1 timevect ydata

objfunction=growthObjective(z,timevect,ydata,flag1,method1);

end
