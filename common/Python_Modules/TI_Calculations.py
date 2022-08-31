
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def processTI(dataTI, restraint, Lsched):
    '''
    Arguments: the TI data, restraint, and lambda schedule
    Function: Calculate the free energy for each lambda value, aggregate the result, and estimate the error
    Returns: The free energies and associated errors as functions of lambda. Both per window and cumulative.
    '''
    dUs = {}
    for key, group in dataTI.groupby('L'):
        dUs[key] = [HW_dUdL(restraint, coord, key) for coord in group.DBC]

    Lsched = np.sort(list(dUs.keys()))
    dL = Lsched[1] - Lsched[0]
    TIperWindow = pd.DataFrame(index=Lsched)
    TIperWindow['dGdL'] = [np.mean(dUs[L])*dL for L in Lsched]
    TIperWindow['error'] = [np.std(dUs[L])*dL for L in Lsched]

    TIcumulative = pd.DataFrame()
    TIcumulative['dG'] = np.cumsum(TIperWindow.dGdL)
    TIcumulative['error'] = np.sqrt(np.divide(np.cumsum(TIperWindow.error**2), np.arange(1,len(TIperWindow)+1)))
    
    return TIperWindow, TIcumulative


def plotTI(cumulative, perWindow, width=8, height=4, color='#0072B2'):
    '''
    Plots TI data with error estimates
    Arguments:
        cumulative: dG data (pandas Series)
        perWindow: dG data (pandas Series)
        width: of figure (inches)
        height: of figure (inches)
        color: of plotted data
    Returns:
        fig, [cumAx, eachAx]: the matplotlib figure and two subplot axes
    '''
    fig, (cumAx,eachAx) = plt.subplots(2,1, sharex='col')

    # Cumulative change in kcal/mol
    cumAx.errorbar(cumulative.index, cumulative.dG, yerr=cumulative.error,marker=None, linewidth=1, color=color, label='Cumulative Change')
    finalEstimate = cumulative.dG[1]
    cumAx.axhline(finalEstimate, linestyle='-', color='gray', label=f'Final Value:\n{np.round(finalEstimate,1)}kcal/mol')
    cumAx.legend()                  
    cumAx.set(ylabel=r'Cumulative $\rm\Delta G_{\lambda}$'+'\n(kcal/mol)')

    # Per-window change in kcal/mol
    eachAx.errorbar(perWindow.index, perWindow.dGdL, marker=None, linewidth=1, yerr=perWindow.error, color=color)
    eachAx.set(ylabel=r'$\rm\Delta G_{\lambda}$'+'\n(kcal/mol)')

    fig.set_figwidth(width)
    fig.set_figheight(height*3)
    fig.tight_layout()
    
    return fig, [cumAx,eachAx] 

def makeHarmonicWall(FC=10, targetFC=0, targetFE=1, upperWalls=1, schedule=None, numSteps=1000, targetEQ=500, name='HW', lowerWalls=None):
    '''
    Makes a "harmonic wall" colvar-like restraint (dict). Includes default values.
    Arguments:
        FC: force constant
        targetFC: target force constant (when lambda=1)
        targetFE: exponent for varying the FC with respect to lambda
        upperWalls: position of the upper wall
        schedule: lambda schedule for a variable FC
        numSteps: the number of steps for each lambda value
        targetEQ: the number of steps discarded to equilibration for each lambda
        name: a name for the harmonic wall (for bookkeeping)
        lowerWalls: the position of the lower wall
    Returns:
        A dict storing all the data
    '''
    HW = {'name':name, 'targetFC':targetFC, 'targetFE':targetFE, 'FC':FC, 'upperWalls':upperWalls, 'schedule':schedule, 'numSteps':numSteps, 'targetEQ':targetEQ, 'lowerWalls':lowerWalls}
    return HW

def HW_U(HW, coord, L):
    '''
    Calculates the potential energy of the harmonic wall.
    Arguments:
        HW: a dict created by makeHarmonicWall
        coord: the value of the colvar
        L: the current lambda value
    Resturns:
        U: the potential energy associated with cv=coord
    '''
    d=0
    if HW['upperWalls'] and coord>HW['upperWalls']:
        d = coord-HW['upperWalls']
    elif HW['lowerWalls'] and coord<HW['lowerWalls']:
        d = coord-HW['lowerWalls']
    
    if d!=0:
        dk = HW['targetFC']-HW['FC']
        la = L**HW['targetFE']
        kL = HW['FC']+la*dk
        U = 0.5*kL*(d**2)
    else:
        U=0
    return U

def HW_dUdL(HW, coord, L):
    '''
    Calculates the potential energy gradient of a harmonic wall.
    Arguments:
        HW: a dict created by makeHarmonicWall
        coord: the value of the colvar
        L: the current lambda value
    Resturns:
        dU: the potential energy gradient associated with cv=coord
    '''
    d=0
    if HW['upperWalls'] and coord>HW['upperWalls']:
        d = coord-HW['upperWalls']
    elif HW['lowerWalls'] and coord<HW['lowerWalls']:
        d = coord-HW['lowerWalls']
    
    if d!=0:
        dk = HW['targetFC']-HW['FC']
        dla = HW['targetFE']*L**(HW['targetFE']-1)
        kL = HW['FC']+dla*dk
        dU = 0.5*kL*(d**2)
    else:
        dU=0
    return dU
