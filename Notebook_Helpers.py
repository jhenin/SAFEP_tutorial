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
    """
    Extracts a number from a file using a regular expression.

    Args:
        regex (str): The regular expression pattern to search for.
        fname (str): The name of the file to search in.
        grp (int, optional): The group number to extract from the regex match. Defaults to 2.

    Returns:
        str: The matched number as a string.

    Example:
        >>> get_num_regex('\d+', 'data.txt')
        '123'
    """
    with open(fname) as fin:
        fstring = fin.read()
        return re.search(regex, fstring).group(grp)
    
# Used in restraint perturbation calculations
uw_regex = '(upperWalls[\ \t]+)(\d+.\d+)'
get_upper_walls = lambda fname: get_num_regex(uw_regex, fname)

def P_bind(K, L):
    """
    Calculate the fraction of binding sites occupied given the dissociation constant and ligand concentration.

    Args:
        K (float): The dissociation constant of the ligand.
        L (float or array-like): The concentration of the ligand.

    Returns:
        float or array-like: The fraction of binding sites occupied, calculated using the formula L / (K + L).
    """
    return L / (K + L)

def Kd(dG, RT):
    """
    Calculate the dissociation constant (Kd) from the Gibbs free energy change (dG) and the product of the gas constant and temperature (RT).

    Args:
        dG (float): The Gibbs free energy change in kcal/mol.
        RT (float): The product of the gas constant (R) and temperature (T) in kcal/mol.

    Returns:
        float: The dissociation constant (Kd) in microMolar (µM).
    """
    return np.exp(dG/RT)*1000000


def plot_titration(ax, concentrations, dG_binding, error_binding, RT):
    """
    Plot the titration curve for a given set of concentrations and binding energies, including a 95% confidence interval.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axis object where the plot will be drawn.
        concentrations (array-like): Array of ligand concentrations in microMolar.
        dG_binding (float): The Gibbs free energy change in kcal/mol for binding.
        error_binding (float): The standard error of the Gibbs free energy change.
        RT (float): The product of the gas constant (R) and temperature (T) in kcal/mol.

    Returns:
        matplotlib.axes.Axes: The axis object with the titration plot.

    """
    K = Kd(dG_binding, RT)

    ax.plot(concentrations, P_bind(K, concentrations), label='Binding Curve')

    P_lower = P_bind(Kd(dG_binding-error_binding*1.96, RT),concentrations)
    P_upper = P_bind(Kd(dG_binding+error_binding*1.96, RT), concentrations)
    ax.fill_between(concentrations, P_lower, P_upper, alpha=0.25, label='95% Confidence Interval')
    plt.xscale('log')
    ax.set_xlabel('Concentration of Phenol '+r'($\mathrm{\mu}$M)', fontsize=20)
    ax.set_ylabel('Fraction of Sites Occupied', fontsize=20)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
    ax.vlines(K, 0, 1, linestyles='dashed', color='black', label='Dissociation Constant')
    ax.legend(loc='lower right', fontsize=20*0.75)

    return ax


@dataclass
class dGData:
    """
    A data class representing dG data for free energy calculations.

    Attributes:
        inpath (Path): The path to the data files.
        filepattern (str): The pattern to match the data file names.
        temperature (float): The temperature in Kelvin.
        name (str): The name of the dG data.
        detectEQ (bool, optional): Whether to detect equilibrium in the data. Defaults to True.
        perWindow (pd.DataFrame, optional): The per-window data. Defaults to an empty DataFrame.
        cumulative (pd.DataFrame, optional): The cumulative data. Defaults to an empty DataFrame.
        dG (float, optional): The calculated dG value. Defaults to 0.
        error (float, optional): The estimated error in the dG value. Defaults to 0.
        RT (float, optional): The product of the gas constant and temperature. Defaults to 0.

    Methods:
        pretty_print_dG(): Returns a formatted string representation of the dG and error values.

    """
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
        """
        Generates a formatted string representing the Gibbs free energy change (ΔG) and its error, using HTML markup for styling.

        Returns:
            str: A string formatted with HTML to display the Gibbs free energy change and its error. The ΔG value is subscripted with the name of the data, and both ΔG and error values are presented in kcal/mol. The string is formatted to display in size 5 font.
        """
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



