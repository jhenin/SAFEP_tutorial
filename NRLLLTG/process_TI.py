import numpy as np
import pandas as pd
import safep
import matplotlib.pyplot as plt
import re
from pymbar.timeseries import detect_equilibration
import argparse

import warnings #Suppress future warnings from pandas.
warnings.simplefilter(action='ignore', category=FutureWarning)

def subsample_traj(dataTI, percent_min, percent_max):
    assert list(dataTI.columns).count('L'), "You must include lambda values for each non-zero energy sample in a column named 'L'"
    
    groups = dataTI.groupby('L')
    trimmed = []
    for key, grp in groups:
        size = len(grp)
        low = int(size*percent_min)
        hi = int(size*percent_max)
        toappend = grp.iloc[low:hi, :]
        trimmed.append(toappend)

    newData = pd.concat(trimmed)
    
    return newData


def read_colvar_traj(colvarsPath, DBC):
    with open(colvarsPath) as f:
        first_line = f.readline()
    columns = re.split(' +', first_line)[1:-1]
    dataTI = pd.read_csv(colvarsPath, delim_whitespace=True, names=columns, comment='#', index_col=0)
    dataTI = dataTI.rename(columns={'E_harmonicwalls2':'DBC_energy'})

    return dataTI


def do_convergence(dataTI, DBC):
    forward = []
    ferr = []
    backward = []
    berr = []
    for x in np.linspace(0.9,0,10):
        subsampled = subsample_traj(dataTI, x, 1)
        TIperWindow, TIcumulative = safep.process_TI(subsampled, DBC, DBC['schedule'])
        backward.append(TIcumulative.dG.iloc[-1])
        berr.append(TIcumulative.error.iloc[-1])
        
        subsampled = subsample_traj(dataTI, 0, 1-x)
        TIperWindow, TIcumulative = safep.process_TI(subsampled, DBC, DBC['schedule'])
        forward.append(TIcumulative.dG.iloc[-1])
        ferr.append(TIcumulative.error.iloc[-1])

    return forward, ferr, backward, berr

def plot_convergence(forward, ferr, backward, berr):
    fig, ax = plt.subplots()
    X = np.linspace(0.1,1,10)
    ax.plot(X, forward, label="forward-time subsampling")
    ax.fill_between(X, np.array(forward)+ferr, np.array(forward)-ferr, alpha = 0.5)
    ax.plot(X, backward, label="backward-time subsampling")
    ax.fill_between(X, np.array(backward)+berr, np.array(backward)-berr, alpha = 0.5)
    ax.set_xlabel("Fraction of non-zero samples")
    ax.set_ylabel("dG (kcal/mol)")
    fig.legend()
    fwd50 = forward[len(forward)//2-1]
    bwd50 = backward[len(backward)//2-1]
    delta50 = fwd50-bwd50
    avg50 = (fwd50+bwd50)/2
    txt = r"$\Delta G_{50} =$"+f"{np.round(delta50,2)}"
    ax.annotate(txt, (0.5, avg50))
    ax.plot([0.5, 0.5], [fwd50, bwd50], color="black")

    return fig, ax

def detect_equilibration_TI(dataTI, DBC):
    groups = dataTI.groupby('L')

    trimmed = []
    for key, grp in groups:
        start, _, eqsamples = detect_equilibration(grp.DBC, nskip=10)
        stride = int(len(grp)/eqsamples)
        toappend = grp.iloc[start::stride, :]
        diff = len(toappend)-eqsamples
        trimmed.append(toappend)

    newData = pd.concat(trimmed)
    newgroups = newData.groupby('L')

    lengths = {}
    for key, grp in newgroups:
        k = (key**DBC['targetFE'])*DBC['targetFC']

        lengths[k] = len(grp)

    return newData, lengths

def plot_samples(lengths):
    fig, ax = plt.subplots()
    toplot = pd.Series(lengths)
    ax.plot(toplot)
    ax.invert_xaxis()
    ax.set_ylabel('Samples')
    ax.set_xlabel(r'$k_\lambda (\frac{kcal}{mol~\AA})$')

    return fig, ax

def setup_TI_analysis(logpath):
    with open(logpath) as logfile:
        wholelog = logfile.read()
        loglines = wholelog.split('\n')
        restraintREGEX = 'harmonicwalls.*initialize'
        restraint = re.search(restraintREGEX, wholelog, re.DOTALL)
        restraintStr = restraint.group(0)
        
        restraint_name = re.search('(.*\")(.*)(\".*\n.*\{\ DBC\ \})', restraintStr).group(2)
        steps_per_lambda = int(re.search('(targetNumSteps\ *=\ *)(\d+)', restraintStr).group(2))
        Lsched = np.float64(re.search('(.*lambdaSchedule\ =\ \{)(.*)(\})', restraintStr).group(2).split(", "))
        FC = float(re.search(f'({restraint_name}.*forceConstant\ *=\ *)(\d+)', restraintStr, re.DOTALL).group(2))
        targetFC = float(re.search('(targetForceConstant\ *=\ *)(\d+)', restraintStr).group(2))
        exponent = float(re.search('(targetForceExponent\ *=\ *)(\d+)', restraintStr).group(2))

        widthRegex = f'({restraint_name}'+'.*upperWalls\ *=\ *{\ )(\d+)(\ })'
        DBCwidth = float(re.search(widthRegex, restraintStr, re.DOTALL).group(2))
        
        eqstepsRegex = f'({restraint_name}'+'.*targetEquilSteps\ *=\ *)(\d+)'
        eqsteps = int(re.search(eqstepsRegex, restraintStr, re.DOTALL).group(2))

    DBC = safep.make_harmonicWall(FC=FC, 
        targetFC=targetFC, 
        targetFE=exponent, 
        upperWalls=DBCwidth, 
        targetEQ=eqsteps, 
        numSteps=steps_per_lambda, 
        name=restraint, 
        schedule=Lsched)

    return DBC, Lsched

def l_from_energy(dataTI):
    guessedLs = [
        safep.guessL(
            U,
            DBC["FC"],
            DBC["targetFC"],
            dbc,
            DBC["upperWalls"],
            DBC["targetFE"],
        )
        for U, dbc in zip(dataTI.DBC_energy, dataTI.DBC)
    ]
    dataTI["L"] = np.round(guessedLs, 6)

    return dataTI

def l_from_sched(dataTI, Lsched):
    Ls = (dataTI.index.values-1)//DBC['numSteps']

    dataLs = np.round([DBC['schedule'][i] for i in Ls], 3)
    dataTI.loc[:,'L'] = dataLs
    dataTI = dataTI.iloc[1:]

    return dataTI

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        prog = "processTI",
        description = "Process the outputs from a NAMD TI calculation."
        )
    parser.add_argument('namdlog')
    parser.add_argument('traj_path')
    parser.add_argument('-o', '--outputDirectory', type=str, default='.')
    parser.add_argument('-d', '--detectEquilibrium', action='store_true')

    args = parser.parse_args()

    print("Reading...")
    DBC, Lsched = setup_TI_analysis(args.namdlog)
    dataTI = read_colvar_traj(args.traj_path, DBC)
    if 'DBC_energy' in dataTI.columns:
        dataTI = l_from_energy(dataTI)
    else:
        dataTI = l_from_sched(dataTI, Lsched)
    
    print("Processing...")
    forward, ferr, backward, berr = do_convergence(dataTI, DBC)
    TIperWindow, TIcumulative = safep.process_TI(dataTI, DBC, Lsched)

    print("Plotting...")
    fig, ax = plot_convergence(forward, ferr, backward, berr)
    plt.savefig(f'{args.outputDirectory}/TI_convergence.pdf')

    fig, ax = safep.plot_TI(TIcumulative, TIperWindow, fontsize=20)
    plt.savefig(f'{args.outputDirectory}/TI_general.pdf')

    if args.detectEquilibrium:
        newData, lengths = detect_equilibration_TI(dataTI, DBC)
        TIperWindow, TIcumulative = safep.process_TI_DBC(newData, DBC)
        fig, ax = plot_samples(lengths)
        plt.savefig(f'{args.outputDirectory}/decorrelated_TI_samples.pdf')

    print("Done")
