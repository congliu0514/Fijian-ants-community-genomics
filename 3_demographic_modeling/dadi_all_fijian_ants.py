from numpy import array
import dadi
import Optimize_Functions
import os
for filename in os.listdir("../input/"):
    sp_name=filename.split(".")[0]
    #pop_size=int(filename.split("-")[1].split(".")[0])
    print (sp_name)
    # Load the data
    data = dadi.Spectrum.from_file("../input/"+filename)


    # Here we are loading the models. 0: Standard neutral model; 1:Exponential change models,
    # and 2:Instantaneous  change model; 3: exponential growth after bottleneck

    func0 = dadi.Demographics1D.snm
    func1 = dadi.Demographics1D.growth
    func2 = dadi.Demographics1D.two_epoch
    #func3 = dadi.Demographics1D.bottlegrowth

    # Optimization parameters for this model.
    # Parameters are: (nu,T). nu: Ratio of contemporary to ancient population size
    #T: Time in the past at which growth began (in units of 2*Na generations)
    pts = [50,70,90]
    upper = [1000, 1000]
    lower = [1e-3, 0]
    p_labels="nu, t"
    reps=[10,20,30,50]
    maxiters=[5,20,125,200]
    folds=[3,2,2,1]

    for i in range(1,3):
        prefix0 = sp_name+"_constant-Replicate_{}".format(i)
        prefix1 = sp_name+"_Exponential-Replicate_{}".format(i)
        prefix2 = sp_name+"_Instantaneous-Replicate_{}".format(i)
        Optimize_Functions.Optimize_Routine(data, pts, prefix0, "constant", func0, 4, 2, fs_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters,folds = folds)
        Optimize_Functions.Optimize_Routine(data, pts, prefix1, "Exponential", func1, 4, 2, fs_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters,folds = folds)
        Optimize_Functions.Optimize_Routine(data, pts, prefix2, "Instantaneous", func2, 4, 2, fs_folded=True, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters,folds = folds)
