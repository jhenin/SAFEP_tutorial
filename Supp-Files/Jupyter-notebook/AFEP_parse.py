#Large datasets can be difficult to parse on a workstation due to inefficiencies in the way data is represented for pymbar. When possible, reduce the size of your dataset.
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from numpy.lib.stride_tricks import sliding_window_view
from scipy.stats import linregress as lr
from scipy.stats import norm
from alchemlyb.preprocessing import subsampling
import os
from scipy.signal import correlate
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit, leastsq


from glob import glob #file regexes
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm #for progress bars
import re #regex
from natsort import natsorted #for sorting "naturally" instead of alphabetically

from alchemlyb.visualisation.dF_state import plot_dF_state

from alchemlyb.parsing import namd
from alchemlyb.estimators import BAR
from alchemlyb.visualisation.dF_state import plot_dF_state
from alchemlyb.visualisation import plot_convergence

import re


def plotGeneral(cumulative, perWindow, RT, width=8, height=4, PDFtype='KDE'):
    fig, ((cumAx, del1),( eachAx, del2),(hystAx, pdfAx)) = plt.subplots(3,2, sharex='col', sharey='row', gridspec_kw={'width_ratios': [2, 1]})

    fig.delaxes(del1)
    fig.delaxes(del2)

    # Cumulative change in kcal/mol
    cumAx.errorbar(cumulative.index, cumulative.BAR.f*RT, yerr=cumulative.BAR.errors, marker=None, linewidth=1)
    cumAx.set(ylabel=r'Cumulative $\rm\Delta G_{\lambda}$'+'\n(kcal/mol)')

    # Per-window change in kcal/mol
    eachAx.errorbar(perWindow.index, perWindow.BAR.df*RT, yerr=perWindow.BAR.ddf, marker=None, linewidth=1)
    eachAx.plot(perWindow.index, perWindow.EXP.dG_f*RT, marker=None, linewidth=1, alpha=0.5)
    eachAx.errorbar(perWindow.index, -perWindow.EXP.dG_b*RT, marker=None, linewidth=1, alpha=0.5)
    eachAx.set(ylabel=r'$\rm\Delta G_{\lambda}$'+'\n(kcal/mol)')

    
    #Hysteresis Plots
    diff = perWindow.EXP['difference']
    hystAx.vlines(perWindow.index, np.zeros(len(perWindow)), diff, label="fwd - bwd", linewidth=2)
    hystAx.set(xlabel=r'$\lambda$', ylabel=r'$\delta_\lambda$ (kcal/mol)', ylim=(-1,1))
    

    
    if PDFtype=='KDE':
        kernel = sp.stats.gaussian_kde(diff)
        pdfX = np.linspace(-1, 1, 1000)
        pdfY = kernel(pdfX)
        pdfAx.plot(pdfY, pdfX, label='KDE')
    elif PDFtype=='PDF':
        pdfX, pdfY = np.histogram(X, nbins=20)
        pdfAx.plot(pdfY, pdfX,  label="Estimated Distribution")
    else:
        raise(f"Error: PDFtype {PDFtype} not recognized")
    
    pdfAx.set(xlabel=PDFtype)

    std = np.std(diff)
    mean = np.average(diff)
    temp = pd.Series(pdfY, index=pdfX)
    mode = temp.idxmax()
    
    textstr = r"$\rm{mode_{pdf}=}$"+f"{np.round(mode,2)}"+"\n"+fr"$\mu$={np.round(mean,2)}"+"\n"+fr"$\sigma$={np.round(std,2)}"
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    pdfAx.text(0.15, 0.95, textstr, transform=pdfAx.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    fig.set_figwidth(width)
    fig.set_figheight(height*3)
    fig.tight_layout()
    
    return fig, [cumAx,eachAx,hystAx,pdfAx] 


def saveUNK(u_nk, filepath):
    u_nk.to_csv(filepath)
    
def readUNK(filepath):
    u_nk = pd.read_csv(filepath)
    u_nk['fep-lambda'] = u_nk['fep-lambda'].astype(str)
    u_nk = u_nk.set_index(['time', 'fep-lambda'])
    
    return u_nk.copy()

#Guess lambda based on file name (last number in the filename divided by 100)
def guessLambda(fname):
    L = int(re.findall(r'\d+', fname)[-1])/100
    return L

#redFEPOUT reads each file in a single pass: keeping track of lambda values and appending each line to an array. 
#The array is cast to a dataframe at the end to avoid appending to a dataframe
def readFEPOUT(fileName, step=1):
    colNames = ["type",'step', 'Elec_l', 'Elec_ldl', 'vdW_l', 'vdW_ldl', 'dE', 'dE_avg', 'Temp', 'dG', 'FromLambda', "ToLambda"]

    data = []

    L = np.nan
    L2 = np.nan
    LIDWS = np.nan
    
    frame = 0
    with open(fileName) as fs:
        for line in fs:
            if line[0] == '#':
                frame = 0
                #print(line)
                Lambda = re.search('LAMBDA SET TO (\d+(\.\d+)*)', line)
                Lambda2 = re.search('LAMBDA2 (\d+(\.\d+)*)', line)
                LambdaIDWS = re.search('LAMBDA_IDWS (\d+(\.\d+)*)', line)
                if Lambda:
                    L = Lambda.group(1)
                    #print(f'L={L}')
                if Lambda2:
                    L2 = Lambda2.group(1)
                    #print(f'L2={L2}')
                if LambdaIDWS:
                    LIDWS = LambdaIDWS.group(1)
                    #print(f'LIDWS={LIDWS}')
            elif frame % step <= 1:
                if np.isnan(L):
                    print("WARNING: lambda is not defined!")
                    L = guessLambda(fileName)
                    print("Guessing lambda to be {L} based on file name.")


                lineList = line.split()
                lineList.append(L)
                if lineList[0] == "FepEnergy:":
                    lineList.append(L2)
                elif lineList[0] == "FepE_back:":
                    lineList.append(LIDWS)
                else:
                    print(f'Unexpected line start: {lineList[0]}')
                    return 0
                data.append(lineList)
                frame = frame + 1
            else:
                frame = frame + 1

            stashL = L
            stashL2 = L2
            stashLIDWS = LIDWS

    fs.close()
    
    df = pd.DataFrame(data).dropna()
    df.columns = colNames
    df = df.iloc[:,1:].astype(float)
    df["window"]=np.mean([df.FromLambda,df.ToLambda], axis=0)
    df["up"]=df.ToLambda>df.FromLambda

    df = df.sort_index()
    return df


# In[5]:


def readFiles(files, step=1):
    fileList = []
    for file in files:
        df = readFEPOUT(file, step)
        fileList.append(df)
    data = pd.concat(fileList)
    
    data.index = data.window
    data["dVdW"] = data.vdW_ldl - data.vdW_l
    
    return data


# In[13]:

def u_nk_fromDF(data, temperature, eqTime, warnings=True):
    from scipy.constants import R, calorie
    beta = 1/(R/(1000*calorie) * temperature) #So that the final result is in kcal/mol
    u_nk = pd.pivot_table(data, index=["step", "FromLambda"], columns="ToLambda", values="dE")
    #u_nk = u_nk.sort_index(level=0).sort_index(axis='columns') #sort the data so it can be interpreted by the BAR estimator
    u_nk = u_nk*beta
    u_nk.index.names=['time', 'fep-lambda']
    u_nk.columns.names = ['']
    u_nk = u_nk.loc[u_nk.index.get_level_values('time')>=eqTime]

    
    #Shift and align values to be consistent with alchemlyb standards
    lambdas = list(set(u_nk.index.get_level_values(1)).union(set(u_nk.columns)))
    lambdas.sort()
    warns = set([])
            
    for L in lambdas:
        try:
            u_nk.loc[(slice(None), L), L] = 0
        except:
            if warnings:
                warns.add(L)
    
    prev = lambdas[0]
    for L in lambdas[1:]:
        try:
            u_nk.loc[(slice(None), L), prev] = u_nk.loc[(slice(None), L), prev].shift(1)
        except:
            if warnings:
                warns.add(L)
            
        prev = L
    
    if len(warns)>0:
        print(f"Warning: lambdas={warns} not found in indices/columns")
    u_nk = u_nk.dropna(thresh=2)
    u_nk = u_nk.sort_index(level=1).sort_index(axis='columns') #sort the data so it can be interpreted by the BAR estimator
    return u_nk


# In[7]:

def readAndProcess(fepoutFiles, temperature, decorrelate, detectEQ):
    from alchemlyb.preprocessing import subsampling

    u_nk = namd.extract_u_nk(fepoutFiles, temperature)
    
    affix=""
    
    if decorrelate:
        print(f"Decorrelating samples. Flag='{decorrelate}'")
        method = 'dE'
        affix = f'{affix}_decorrelated_{method}'
        groups = u_nk.groupby('fep-lambda')
        decorr = pd.DataFrame([])
        for key, group in groups:
            test = subsampling.decorrelate_u_nk(group, method)
            decorr = pd.concat([decorr, test])
        u_nk = decorr
    else:
        affix = f'{affix}_unprocessed'
    
    if detectEQ:
        print("Detecting Equilibrium")
        affix = f"{affix}_AutoEquilibrium"
        groups = u_nk.groupby('fep-lambda')
        EQ = pd.DataFrame([])
        for key, group in groups:
            group = group[~group.index.duplicated(keep='first')]
            test = subsampling.equilibrium_detection(group, group.dropna(axis=1).iloc[:,-1])
            EQ = pd.concat([EQ, test])
        u_nk = EQ
    else:
        affix=f"{affix}_HardEquilibrium"

    return u_nk, affix


def get_dG(u_nk):
    #the data frame is organized from index level 1 (fep-lambda) TO column
    #dG will be FROM column TO index
    groups = u_nk.groupby(level=1)
    dG=pd.DataFrame([]) 
    for name, group in groups:
        dG[name] = np.log(np.mean(np.exp(-1*group)))
        dG = dG.copy() # this is actually faster than having a fragmented dataframe
        
    return dG


def doEstimation(u_nk, method='both'):
    u_nk = u_nk.sort_index(level=1)
    cumulative = pd.DataFrame()
    perWindow = pd.DataFrame()
    if method=='both' or method=='BAR':
        bar = BAR()
        bar.fit(u_nk)
        ls, l_mids, fs, dfs, ddfs, errors = get_BAR(bar)
        
        cumulative[('BAR', 'f')] = fs
        cumulative[('BAR', 'errors')] = errors
        cumulative.index = ls

        perWindow[('BAR','df')] = dfs
        perWindow[('BAR', 'ddf')] = ddfs
        perWindow.index = l_mids
        
    if method=='both' or method=='EXP':
        expl, expmid, dG_fs, dG_bs = get_EXP(u_nk)

        cumulative[('EXP', 'ff')] = np.insert(np.cumsum(dG_fs),0,0)
        cumulative[('EXP', 'fb')] = np.insert(-np.cumsum(dG_bs),0,0)
        cumulative.index = expl 
        
        perWindow[('EXP','dG_f')] = dG_fs
        perWindow[('EXP','dG_b')] = dG_bs
        perWindow[('EXP', 'difference')] = np.array(dG_fs)+np.array(dG_bs)        
        perWindow.index = expmid
        
    
    perWindow.columns = pd.MultiIndex.from_tuples(perWindow.columns)
    cumulative.columns = pd.MultiIndex.from_tuples(cumulative.columns)
    
    return perWindow.copy(), cumulative.copy()




def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w


# Subsamples a u_nk dataframe using percentiles [0-100] of data instead of absolute percents
def subSample(unkGrps, lowPct, hiPct):
    partial = []
    for key, group in unkGrps:
        idcs = group.index.get_level_values(0)
        
        lowBnd = np.percentile(idcs, lowPct, method='closest_observation')
        hiBnd = np.percentile(idcs, hiPct, method='closest_observation')
        mask = np.logical_and(idcs<=hiBnd, idcs>=lowBnd) 
        sample = group.loc[mask]
        if len(sample)==0:
            print(f"ERROR: no samples in window {key}")
            print(f"Upper bound: {hiBnd}\nLower bound: {lowBnd}")
            raise
            
        partial.append(sample)

    partial = pd.concat(partial)
    
    return partial


# altConvergence splits the data into percentile blocks. Inspired by block averaging
def altConvergence(u_nk, nbins):
    groups = u_nk.groupby('fep-lambda')

    #return data_list
    
    forward = []
    forward_error = []
    backward = []
    backward_error = []
    num_points = nbins
    for i in range(1, num_points+1):
        # forward
        partial = subSample(groups, 100*(i-1)/num_points, 100*i/num_points)
        estimate = BAR().fit(partial)
        l, l_mid, f, df, ddf, errors = get_BAR(estimate)
        
        forward.append(f.iloc[-1])
        forward_error.append(errors[-1])

    return np.array(forward), np.array(forward_error)

def doConvergence(u_nk, tau=1, num_points=10):
    groups = u_nk.groupby('fep-lambda')

    #return data_list
    
    forward = []
    forward_error = []
    backward = []
    backward_error = []
    for i in range(1, num_points+1):
        # forward
        partial = subSample(groups, 0, 100*i/num_points)
        estimate = BAR().fit(partial)
        l, l_mid, f, df, ddf, errors = get_BAR(estimate)
        
        forward.append(f.iloc[-1])
        forward_error.append(errors[-1])
        
        partial = subSample(groups, 100*(1-i/num_points), 100)
        estimate = BAR().fit(partial)
        l, l_mid, f, df, ddf, errors = get_BAR(estimate)
        
        backward.append(f.iloc[-1])
        backward_error.append(errors[-1])

    return np.array(forward), np.array(forward_error), np.array(backward), np.array(backward_error)



def doConvPlot(ax, X, fs, ferr, fwdColor, label=None):
    ax.errorbar(X, fs, yerr=ferr, marker=None, linewidth=1, color=fwdColor, markerfacecolor='white', markeredgewidth=1, markeredgecolor=fwdColor, ms=5, label=label)
    return ax



def convergencePlot(theax, fs, ferr, bs, berr, fwdColor='#0072B2', bwdColor='#D55E00', lgndF=None, lgndB=None):
    if not lgndF:
        lgndF=fwdColor
        lgndB=bwdColor
        
        
    lower = fs[-1]-ferr[-1]
    upper = fs[-1]+ferr[-1]
    theax.fill_between([0,1],[lower, lower], [upper, upper], color=bwdColor, alpha=0.25)
    theax.errorbar(np.arange(len(fs))/len(fs)+0.1, fs, yerr=ferr, marker='o', linewidth=1, color=fwdColor, markerfacecolor='white', markeredgewidth=1, markeredgecolor=fwdColor, ms=5)
    theax.errorbar(np.arange(len(bs))/len(fs)+0.1, bs, yerr=berr, marker='o', linewidth=1, color=bwdColor, markerfacecolor='white', markeredgewidth=1, markeredgecolor=bwdColor, ms=5, linestyle='--')


    
    theax.xaxis.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    
    finalMean = fs[-1]
    theax.axhline(y= finalMean, linestyle='-.', color='gray')
    theax.plot(0, finalMean, linewidth=1, color=lgndF, label='Forward Time Sampling')
    theax.plot(0, finalMean, linewidth=1, color=lgndB, linestyle='--', label='Backward Time Sampling')
    
    return theax















#Cannonical convergence plot
def convergence_plot(u_nk, tau=1, units='kT', RT=0.59):
    forward, forward_error, backward, backward_error = doConvergence(u_nk, num_points=10)

    if units=='kcal/mol':
        forward = forward*RT
        forward_error = forward_error*RT
        backward = backward*RT
        backward_error = backward_error*RT
    
    ax = plot_convergence(forward, forward_error, backward, backward_error)
    
    if units=='kcal/mol':
        ax.set(ylabel=r'$\rm\Delta G$'+'\n(kcal/mol)')

    return plt.gca()

def get_BAR(bar):
    
    # Extract data for plotting
    states = bar.states_

    f = bar.delta_f_.iloc[0,:] # dataframe
    l = np.array([float(s) for s in states])
    # lambda midpoints for each window
    l_mid = 0.5*(l[1:] + l[:-1])

    # FE differences are off diagonal
    df = np.diag(bar.delta_f_, k=1)
    

    # error estimates are off diagonal
    ddf = np.array([bar.d_delta_f_.iloc[i, i+1] for i in range(len(states)-1)])

    # Accumulate errors as sum of squares
    errors = np.array([np.sqrt((ddf[:i]**2).sum()) for i in range(len(states))])
    
    
    return l, l_mid, f, df, ddf, errors

def get_EXP(u_nk):
    #the data frame is organized from index level 1 (fep-lambda) TO column
    #dG will be FROM column TO index
    groups = u_nk.groupby(level=1)
    dG=pd.DataFrame([])
    for name, group in groups:
        dG[name] = np.log(np.mean(np.exp(-1*group)))

    dG_f=np.diag(dG, k=1)
    dG_b=np.diag(dG, k=-1)

    l=dG.columns.to_list()
    l_mid = np.mean([l[1:],l[:-1]], axis=0)

    return l, l_mid, dG_f, dG_b

def fb_discrepancy_plot(l_mid, dG_f, dG_b):
    plt.vlines(l_mid, np.zeros(len(l_mid)), dG_f + np.array(dG_b), label="fwd - bwd", linewidth=3)
    plt.legend()
    plt.title('Fwd-bwd discrepancies by lambda')
    plt.xlabel('Lambda')
    plt.ylabel('Diff. in delta-G')
    #Return the figure
    return plt.gca() 

def fb_discrepancy_hist(dG_f, dG_b):
    plt.hist(dG_f + np.array(dG_b));
    plt.title('Distribution of fwd-bwd discrepancies')
    plt.xlabel('Difference in delta-G')
    plt.ylabel('Count')
    #return the figure
    return plt.gca()


#Light-weight exponential estimator
def get_dG_fromData(data, temperature):
    from scipy.constants import R, calorie
    beta = 1/(R/(1000*calorie) * temperature) #So that the final result is in kcal/mol
    
    groups = data.groupby(level=0)
    dG=[]
    for name, group in groups:
        isUp = group.up
        dE = group.dE
        toAppend = [name, -1*np.log(np.mean(np.exp(-beta*dE[isUp]))), 1]
        dG.append(toAppend)
        toAppend=[name, -1*np.log(np.mean(np.exp(-beta*dE[~isUp]))), 0]
        dG.append(toAppend)
    
    dG = pd.DataFrame(dG, columns=["window", "dG", "up"])
    dG = dG.set_index('window')
    
    dG_f = dG.loc[dG.up==1] 
    dG_b = dG.loc[dG.up==0]

    dG_f = dG_f.dG.dropna()
    dG_b = dG_b.dG.dropna()

    return dG_f, dG_b


# Additional utility functions:

#Functions for bootstrapping estimates and generating confidence intervals
def bootStrapEstimate(u_nk, estimator='BAR', iterations=100, schedule=[10,20,30,40,50,60,70,80,90,100]):
    groups = u_nk.groupby('fep-lambda')

    if estimator == 'EXP':
        dGfs = {}
        dGbs = {}
        alldGs = {}
    elif estimator == 'BAR':
        dGs = {}
        errs = {}
    else:
        raise ValueError(f"unknown estimator: {estimator}")

    for p in schedule:
        Fs = []
        Bs = []
        fs = []
        Gs = []
        #rs = []
        for i in np.arange(iterations):
            sampled = pd.DataFrame([])
            for key, group in groups:
                N = int(p*len(group)/100)
                if N < 1:
                    N=1
                rows = np.random.choice(len(group), size=N)
                test = group.iloc[rows,:]
                sampled = pd.concat([sampled, test])
            if estimator == 'EXP':
                l, l_mid, dG_f, dG_b = get_EXP(pd.DataFrame(sampled))
                F = np.sum(dG_f)
                B = np.sum(-dG_b)
                Fs.append(F)
                Bs.append(B)
                Gs.append(np.mean([F,B]))
            elif estimator == 'BAR':
                tmpBar = BAR()
                tmpBar.fit(sampled)
                l, l_mid, f, df, ddf, errors = get_BAR(tmpBar)
                fs.append(f.values[-1])
                #rs.append(errors[-1])

        if estimator == 'EXP':
            dGfs[p] = Fs
            dGbs[p] = Bs
            alldGs[p] = Gs
        else:
            dGs[p] = fs
            #errs[p] = rs

    if estimator == 'EXP':
        fwd = pd.DataFrame(dGfs).melt().copy()
        bwd = pd.DataFrame(dGbs).melt().copy()
        alldGs = pd.DataFrame(alldGs).melt().copy()
        return (alldGs, fwd, bwd)
    else:
        alldGs = pd.DataFrame(dGs).melt().copy()
        #allErrors = pd.DataFrame(errs).melt().copy()
        return alldGs
    
def getLimits(allSamples):
    groups = allSamples.groupby('variable')
    means = []
    errors = []
    for key, group in groups:
        means.append(np.mean(group.value))
        errors.append(np.std(group.value))

    upper = np.sum([[x*1 for x in errors],means], axis=0)
    lower = np.sum([[x*(-1) for x in errors],means], axis=0)
    
    return (upper, lower, means)

def getEmpiricalCI(allSamples, CI=0.95):
    groups = allSamples.groupby('variable')

    uppers=[]
    lowers=[]
    means=[]
    for key, group in groups:
        uppers.append(np.sort(group.value)[round(len(group)*CI)])
        lowers.append(np.sort(group.value)[round(len(group)*(1-CI))])
        means.append(np.mean(group.value))

    return (uppers, lowers, means)


# Estimate the probability density distribution from the moving slope of a CDF. i.e. using the values X and their cumulative density FX
def getMovingAveSlope(X,FX,window):
    slopes = []
    Xwindowed = sliding_window_view(X, window)
    FXwindowed = sliding_window_view(FX, window)
   
    for i in np.arange(len(Xwindowed)):
        Xwindow = Xwindowed[i]
        FXwindow = FXwindowed[i]
        result = lr(Xwindow, FXwindow)
        m = result.slope
        slopes.append(m)
    return slopes

# Calculate the coefficient of determination:
def GetRsq(X, Y, Yexpected):
    residuals = Y-Yexpected
    SSres = np.sum(residuals**2)
    SStot = np.sum((X-np.mean(X))**2)
    R = 1-SSres/SStot
    R


from scipy.special import erfc
from scipy.optimize import curve_fit as scipyFit
from scipy.stats import skew
#Wrapper for fitting the normal CDF
def cumFn(x, m, s):
    r = norm.cdf(x, m, s)
    return r

def pdfFn(x,m,s):
    r = norm.pdf(x,m,s)
    return r

#Calculate the PDF of the discrepancies
def getPDF(dG_f, dG_b, DiscrepancyFitting='LS', dx=0.01, binNum=20):
    diff = dG_f + np.array(dG_b)
    diff.sort()
    X = diff
    Y = np.arange(len(X))/len(X)

    #fit a normal distribution to the existing data
    if DiscrepancyFitting == 'LS':
        fitted = scipyFit(cumFn, X, Y)[0] #Fit norm.cdf to (X,Y)
    elif DiscrepancyFitting == 'ML':
        fitted = norm.fit(X) # fit a normal distribution to X
    else:
        raise("Error: Discrepancy fitting code not known. Acceptable values: ML (maximum likelihood) or LS (least squares)")
    discrepancies = dG_f + np.array(dG_b)

    pdfY, pdfX = np.histogram(discrepancies, bins=binNum, density=True)
    pdfX = (pdfX[1:]+pdfX[:-1])/2

    pdfXnorm  = np.arange(np.min(X), np.max(X), dx)
    pdfYnorm = norm.pdf(pdfXnorm, fitted[0], fitted[1])

    pdfYexpected = norm.pdf(pdfX, fitted[0], fitted[1])
           
    return X, Y, pdfX, pdfY, fitted, pdfXnorm, pdfYnorm, pdfYexpected


if __name__ == '__main__':

    import argparse
    import os
  
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--path', type=str, help='The absolute path to the directory containing the fepout files', default='.')
    parser.add_argument('--fepoutre', type=str, help='A regular expression that matches the fepout files of interest.', default='*.fep*')
    parser.add_argument('--temperature', type=float, help='The temperature at which the FEP was run.')
    parser.add_argument('--decorrelate', type=bool, help='Flag to determine whether or not to decorrelate the data. (1=decorrelate, 0=use all data)', default=0)
    parser.add_argument('--detectEQ', type=bool, help='Flag for automated equilibrium detection.', default=0)
    parser.add_argument('--fittingMethod', type=str, help='Method for fitting the forward-backward discrepancies (hysteresis). LS=least squares, ML=maximum likelihood Default: LS', default='LS')
    parser.add_argument('--maxSize', type=float, help='Maximum total file size in GB. This is MUCH less than the required RAM. Default: 1', default=1)
    parser.add_argument('--makeFigures', type=bool, help='Run additional diagnostics and save figures to the directory. default: False', default=0)

    args = parser.parse_args()

    path = args.path
    filename = args.fepoutre
    maxSize = args.maxSize
    temperature = args.temperature
    decorrelate = (args.decorrelate==1)
    detectEQ = (args.detectEQ==1)
    DiscrepancyFitting = args.fittingMethod
    
    RT = 0.00198720650096 * temperature # ca. 0.59kcal/mol


    fepoutFiles = glob(path+filename)
    print(f"Processing: {path+filename}")

    totalSize = 0
    for file in fepoutFiles:
        totalSize += os.path.getsize(file)
    print(f"Fepout Files: {len(fepoutFiles)}\n")

    if (totalSize/10**9)>maxSize:
        print(f"Error: desired targets (Total size:{np.round(totalSize/10**9, 2)}GB) exceed your max size. Either increase your maximum acceptable file size, or use the 'Extended' notebook")
        raise


    print(f'DetectEQ: {detectEQ}')
    print(f'Decorr: {decorrelate}')
    u_nk, affix = readAndProcess(fepoutFiles, temperature, decorrelate, detectEQ)


    u_nk = u_nk.sort_index(level=1)
    bar = BAR()
    bar.fit(u_nk)
    l, l_mid, f, df, ddf, errors = get_BAR(bar)
    changeAndError = f'\u0394G = {np.round(f.iloc[-1]*RT, 1)}\u00B1{np.round(errors[-1], 3)} kcal/mol'
    print(changeAndError)

    if args.makeFigures == 1:
        # Cumulative change in kT
        plt.errorbar(l, f, yerr=errors, marker='.')
        plt.xlabel('lambda')
        plt.ylabel('DeltaG(lambda) (kT)')
        plt.title(f'Cumulative dG with accumulated errors {affix}\n{changeAndError}')
        plt.savefig(f'{path}dG_cumulative_kT_{affix}.png', dpi=600)
        plt.clf()

        # Per-window change in kT
        plt.errorbar(l_mid, df, yerr=ddf, marker='.')
        plt.xlabel('lambda')
        plt.ylabel('Delta G per window (kT)')
        plt.title(f'Per-Window dG with individual errors {affix}')
        plt.savefig(f'{path}dG_{affix}.png', dpi=600)
        plt.clf()

        # Per-window change in kT
        plt.errorbar(l[1:-1], np.diff(df), marker='.')
        plt.xlabel('lambda (L)')
        plt.ylabel("dG'(L)")
        plt.title(f'derivative of dG {affix}')
        plt.savefig(f'{path}dG_prime_{affix}.png', dpi=600)
        plt.clf()

        ####
        try:
            convergence_plot(u_nk, l)
            plt.title(f'Convergence {affix}')
            plt.savefig(f'{path}convergence_{affix}.png', dpi=600)
            plt.clf()
        except:
            print("Failed to generate convergence plot. Probably due to too few samples after decorrelation.")

        ####
        l, l_mid, dG_f, dG_b = get_EXP(u_nk)
        plt.vlines(l_mid, np.zeros(len(l_mid)), dG_f + np.array(dG_b), label="fwd - bwd", linewidth=2)

        plt.legend()
        plt.title(f'Fwd-bwd discrepancies by lambda {affix}')
        plt.xlabel('Lambda')
        plt.ylabel('Diff. in delta-G')
        plt.savefig(f'{path}discrepancies_{affix}.png', dpi=600)
        plt.clf()


        ###
        #Do residual fitting
        ###
        X, Y, pdfX, pdfY, fitted, pdfXnorm, pdfYnorm, pdfYexpected = getPDF(dG_f, dG_b)
        
        #plot the data
        fig, (pdfAx, pdfResid) = plt.subplots(2, 1, sharex=True)
        plt.xlabel('Difference in delta-G')
        
        pdfAx.plot(pdfX, pdfY,  label="Estimated Distribution")
        pdfAx.set_ylabel("PDF")
        pdfAx.plot(pdfXnorm, pdfYnorm, label="Fitted Normal Distribution", color="orange")

        #pdf residuals
        pdfResiduals = pdfY-pdfYexpected
        pdfResid.plot(pdfX, pdfResiduals)
        pdfResid.set_ylabel("PDF residuals") 

        fig.set_figheight(10)
        if DiscrepancyFitting == 'LS':
            pdfAx.title.set_text(f"Least squares fitting of cdf(fwd-bkwd)\nSkewness: {np.round(skew(X),2)}\nFitted parameters: Mean={np.round(fitted[0],3)}, Stdv={np.round(fitted[1],3)}\nPopulation parameters: Mean={np.round(np.average(X),3)}, Stdv={np.round(np.std(X),3)}")
            plt.savefig(f"{path}LeastSquares_pdf_{affix}.png", dpi=600)
        elif DiscrepancyFitting == 'ML':
            pdfAx.title.set_text(f"Maximum likelihood fitting of fwd-bkwd\nFitted parameters: Mean={np.round(fitted[0],3)}, Stdv={np.round(fitted[1],3)}\nPopulation parameters: Mean={np.round(np.average(X),3)}, Stdv={np.round(np.std(X),3)}")
            plt.savefig(f"{path}MaximumLikelihood_pdf_{affix}.png", dpi=600)
        plt.clf()

