"""
Helper functions to support the SAFEP tutorial notebook.
Some of these will eventually end up in the SAFEP package.
"""
import numpy as np
import pandas as pd
import scipy as sp
from dataclasses import dataclass
from pathlib import Path
import safep
import re


@dataclass
class FEPData:
    inpath: Path
    filepattern: str
    temperature: float
    detectEQ: bool
    name: str
    perWindow: pd.DataFrame = pd.DataFrame(None)
    cumulative: pd.DataFrame = pd.DataFrame(None)
    forward: pd.Series = pd.Series(None)
    forward_error: pd.Series = pd.Series(None)
    backward: pd.Series = pd.Series(None)
    backward_error: pd.Series = pd.Series(None)
    dG: float = 0
    error: float = 0
    RT: float = 0

    def process(self):
        R = sp.constants.R/(1000*sp.constants.calorie) # gas constant in kcal/(mol K)
        self.RT = R * self.temperature # RT in kcal/mol

        fepoutFiles = self.inpath.glob(self.filepattern) #Resolve any naming regular expressions
        u_nk_site = safep.read_and_process(fepoutFiles, self.temperature, decorrelate=False, detectEQ=self.detectEQ) #u_nk stores the fep data
        self.perWindow, self.cumulative = safep.do_estimation(u_nk_site) #Run the BAR estimator on the fep data
        self.forward, self.forward_error, self.backward, self.backward_error = safep.do_convergence(u_nk_site) #Used later in the convergence plot

        self.dG = np.round(self.cumulative.BAR.f.iloc[-1]*self.RT, 1)
        self.error = np.round(self.cumulative.BAR.errors.iloc[-1]*self.RT, 1)

    def pretty_print_dG(self):
        #\u0394 == Delta
        change_mkd_site = f'\u0394G<sub>{self.name}</sub> = {self.dG} kcal/mol'
        error_mkd_site = f'PyMBAR estimated error: {self.error} kcal/mol'
        mkd_string = '<font size=5>{}</font><br/><font size=5>{}</font><br/>'.format(change_mkd_site, error_mkd_site)
        return mkd_string


def get_num_regex(regex, fname, grp=2):
    with open(fname) as fin:
        fstring = fin.read()
        return re.search(regex, fstring).group(grp)
    
# Used in restraint perturbation calculations
uw_regex = '(upperWalls[\ \t]+)(\d+.\d+)'
get_upper_walls = lambda fname: get_num_regex(uw_regex, fname)




