"""
Helper functions to support the SAFEP tutorial notebook.
Some of these will eventually end up in the SAFEP package.
"""
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from dataclasses import dataclass
from pathlib import Path
import safep
import re


def get_num_regex(regex, fname, grp=2):
    with open(fname) as fin:
        fstring = fin.read()
        return re.search(regex, fstring).group(grp)
    
# Used in restraint perturbation calculations
uw_regex = '(upperWalls[\ \t]+)(\d+.\d+)'
get_upper_walls = lambda fname: get_num_regex(uw_regex, fname)

def P_bind(K, L):
    return L/(K+L)

def Kd(dG, RT):
    return np.exp(dG/RT)*1000000

@dataclass
class dGData:
    inpath: Path
    filepattern: str
    temperature: float
    name: str
    detectEQ: bool = True
    perWindow: pd.DataFrame = pd.DataFrame(None)
    cumulative: pd.DataFrame = pd.DataFrame(None)
    dG: float = 0
    error: float = 0
    RT: float = 0

    def pretty_print_dG(self):
        #\u0394 == Delta
        change_mkd_site = f'\u0394G<sub>{self.name}</sub> = {self.dG} kcal/mol'
        error_mkd_site = f'PyMBAR estimated error: {self.error} kcal/mol'
        mkd_string = '<font size=5>{}</font><br/><font size=5>{}</font><br/>'.format(change_mkd_site, error_mkd_site)
        return mkd_string


@dataclass
class TIData(dGData):
    eqtime: int = 1000
    upper_walls: float=None
    nLambdas: int=41
    Lsched: list = None
    harmonic_wall: dict = None
    force_constant: float=0
    target_force_constant: float=200
    force_exponent: float=6
    num_steps: int = 300000
    data: pd.DataFrame = pd.DataFrame(None)

    def read(self):
        # Setup and processing of colvars data
        infile = list(self.inpath.glob(self.filepattern))[0]
        with open(infile) as f:
            first_line = f.readline()
        columns = re.split(' +', first_line)[1:-1]
        self.data = pd.read_csv(infile,
                             delim_whitespace=True,
                             names=columns,
                             comment='#',
                             index_col=0)
        self.data = self.data[self.data.index>=self.eqtime][1:]
        self.data.index = self.data.index-self.eqtime

        return

    def process(self):
        if not self.Lsched:
            self.Lsched = np.linspace(1,0,self.nLambdas)

        self.harmonic_wall = safep.make_harmonicWall(FC=self.force_constant, 
                                      targetFC=self.target_force_constant,
                                      targetFE=self.force_exponent,
                                      upperWalls=self.upper_walls,
                                      targetEQ=self.eqtime,
                                      numSteps=self.num_steps,
                                      name=self.name,
                                      schedule=self.Lsched)
        
        Ls = (self.data.index.values-1)//self.harmonic_wall['numSteps']
        Ls[0] = 0

        # This is a small hack in case there are extra samples for the last window
        Ls[Ls==self.nLambdas] = self.nLambdas-1 

        dataLs = np.round([self.harmonic_wall['schedule'][i] for i in Ls], 3)
        self.data.loc[:,'L'] = dataLs
        self.data = self.data.iloc[1:]

        self.perWindow, self.cumulative = safep.process_TI(self.data,
                                                           self.harmonic_wall,
                                                           self.Lsched)
        self.dG = np.round(self.cumulative['dG'][1], 1)
        self.error = np.round(self.cumulative['error'][1], 1)

        return

@dataclass
class FEPData(dGData):
    forward: pd.Series = pd.Series(None)
    forward_error: pd.Series = pd.Series(None)
    backward: pd.Series = pd.Series(None)
    backward_error: pd.Series = pd.Series(None)

    def process(self):
        R = sp.constants.R/(1000*sp.constants.calorie) # gas constant in kcal/(mol K)
        self.RT = R * self.temperature # RT in kcal/mol

        fepoutFiles = self.inpath.glob(self.filepattern) #Resolve any naming regular expressions
        u_nk_site = safep.read_and_process(fepoutFiles, self.temperature, decorrelate=False, detectEQ=self.detectEQ) #u_nk stores the fep data
        self.perWindow, self.cumulative = safep.do_estimation(u_nk_site) #Run the BAR estimator on the fep data
        self.forward, self.forward_error, self.backward, self.backward_error = safep.do_convergence(u_nk_site) #Used later in the convergence plot

        self.dG = np.round(self.cumulative.BAR.f.iloc[-1]*self.RT, 1)
        self.error = np.round(self.cumulative.BAR.errors.iloc[-1]*self.RT, 1)

        return

    def general_plot(self, width, height, cumulative_ylim, perwindow_ylim):
        fig, axes = safep.plot_general(
                    self.cumulative,
                    cumulative_ylim, 
                    self.perWindow, 
                    perwindow_ylim, 
                    self.RT, 
                    width=width, 
                    height=height, 
                    pdf_type='KDE', 
                    fontsize=20)
        return fig, axes

    def convergence_plot(self, width, height, fontsize=20):
        fig, convAx = plt.subplots(1,1)
        convAx = safep.convergence_plot(convAx,
                                        self.forward*self.RT, 
                                        self.forward_error*self.RT,
                                        self.backward*self.RT,
                                        self.backward_error*self.RT, 
                                        fontsize=fontsize)

        fig.set_figwidth(width)
        fig.set_figheight(height)

        return fig, convAx



