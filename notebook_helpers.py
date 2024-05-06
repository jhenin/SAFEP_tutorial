"""
Helper functions to support the SAFEP tutorial notebook.
Some of these will eventually end up in the SAFEP package.

Linting comments:
These are standard abbreviations:
dG = delta Gibb's free energy
RT = gas constant * temperature
"""

import re
from pathlib import Path
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import safep


def get_num_regex(regex, fname, grp=2):
    """
    Extracts a number from a file using a regular expression.

    Args:
        regex (str): The regular expression pattern to search for.
        fname (str): The name of the file to search in.
        grp (int, optional): The group number to extract from the regex match. Defaults to 2.

    Returns:
        str: The matched number as a string.
    """
    with open(fname, "r", encoding="UTF8") as fin:
        fstring = fin.read()
        found = re.search(regex, fstring)
        if found is not None:
            toreturn = found.group(grp)
        else:
            toreturn = None
        return toreturn


# Used in restraint perturbation calculations
UW_REGEX = r"(upperWalls[\ \t]+)(\d+.\d+)"


def get_upper_walls(fname):
    """
    get the upper walls from a colvars config file
    """
    return get_num_regex(UW_REGEX, fname)


def binding_probability(dissociation_constant, concentrations):
    """
    Calculate the fraction of binding sites occupied given the dissociation constant
    and ligand concentration.

    Args:
        dissociation_constant (float): The dissociation constant of the ligand.
        concentrations (float or array-like): The concentration of the ligand.

    Returns:
        float or array-like: The fraction of binding sites occupied, calculated
        using the formula L / (K + L).
    """
    return concentrations / (dissociation_constant + concentrations)


def get_dissociation_constant(delta_g, RT):
    """
    Calculate the dissociation constant (Kd) from the Gibbs free energy change (delta_g)
    and the product of the gas constant and temperature (RT).

    Args:
        delta_g (float): The Gibbs free energy change in kcal/mol.
        RT (float): The product of the gas constant (R) and temperature (T) in kcal/mol.

    Returns:
        float: The dissociation constant (Kd) in microMolar (µM).
    """
    return np.exp(delta_g / RT) * 1000000


def plot_titration(axis, concentrations, binding_free_energy, error_binding, RT):
    """
    Plot the titration curve for a given set of concentrations and binding energies,
    including a 95% confidence interval.

    Args:
        axis (matplotlib.axes.Axes): The matplotlib axis object where the plot will be drawn.
        concentrations (array-like): Array of ligand concentrations in microMolar.
        delta_g_binding (float): The Gibbs free energy change in kcal/mol for binding.
        error_binding (float): The standard error of the Gibbs free energy change.
        RT (float): The product of the gas constant (R) and temperature (T) in kcal/mol.

    Returns:
        matplotlib.axes.Axes: The axis object with the titration plot.

    """
    k_dissociation = get_dissociation_constant(binding_free_energy, RT)

    axis.plot(
        concentrations,
        binding_probability(k_dissociation, concentrations),
        label="Binding Curve",
    )

    probability_lower_bound = binding_probability(
        get_dissociation_constant(binding_free_energy - error_binding * 1.96, RT),
        concentrations,
    )
    probability_upper_bound = binding_probability(
        get_dissociation_constant(binding_free_energy + error_binding * 1.96, RT),
        concentrations,
    )
    axis.fill_between(
        concentrations,
        probability_lower_bound,
        probability_upper_bound,
        alpha=0.25,
        label="95% Confidence Interval",
    )
    plt.xscale("log")
    axis.set_xlabel("Concentration of Phenol " + r"($\mathrm{\mu}$M)", fontsize=20)
    axis.set_ylabel("Fraction of Sites Occupied", fontsize=20)
    axis.set_xticklabels(axis.get_xticklabels(), fontsize=16)
    axis.set_yticklabels(axis.get_yticklabels(), fontsize=16)
    axis.vlines(
        k_dissociation,
        0,
        1,
        linestyles="dashed",
        color="black",
        label="Dissociation Constant",
    )
    axis.legend(loc="lower right", fontsize=20 * 0.75)

    return axis


@dataclass
class DeltaGData:
    """
    A data class representing dG data for free energy calculations.

    Attributes:
        inpath (Path): The path to the data files.
        filepattern (str): The pattern to match the data file names.
        temperature (float): The temperature in Kelvin.
        name (str): The name of the dG data.
        detect_equilibrium (bool, optional): Whether to detect equilibrium in the data.
            Defaults to True.
        per_window (pd.DataFrame, optional): The per-window data.
            Defaults to an empty DataFrame.
        cumulative (pd.DataFrame, optional): The cumulative data.
            Defaults to an empty DataFrame.
        dG (float, optional): The calculated dG value. Defaults to 0.
        error (float, optional): The estimated error in the dG value. Defaults to 0.
        RT (float, optional): The product of the gas constant and temperature.
            Defaults to 0.

    Methods:
        pretty_print_delta_g(): Returns a formatted string representation of the dG
            and error values.

    """

    # pylint: disable=too-many-instance-attributes
    # we might deepen the class heirarchy later, but it works for now.

    inpath: Path
    filepattern: str
    temperature: float
    name: str
    detect_equilibrium: bool = True
    per_window: pd.DataFrame = field(default_factory=pd.DataFrame)
    cumulative: pd.DataFrame = field(default_factory=pd.DataFrame)
    delta_g: float = 0
    error: float = 0
    RT: float = 0

    def pretty_print_delta_g(self):
        """
        Generates a formatted string representing the Gibbs free energy change
        (ΔG) and its error, using HTML markup for styling.

        Returns:
            str: A string formatted with HTML to display the Gibbs free energy
            change and its error. The ΔG value is subscripted with the name of the data,
            and both ΔG and error values are presented in kcal/mol. The string is
            formatted to display in size 5 font.
        """
        # \u0394 == Delta
        change_mkd_site = f"\u0394G<sub>{self.name}</sub> = {self.delta_g} kcal/mol"
        error_mkd_site = f"PyMBAR estimated error: {self.error} kcal/mol"
        mkd_string = f"<font size=5>{change_mkd_site}</font><br/>"
        mkd_string += f"<font size=5>{error_mkd_site}</font><br/>"
        return mkd_string


@dataclass
class TIData(DeltaGData):
    """
    A data class representing TIData for free energy calculations.

    Attributes:
        inpath (Path): The path to the data files.
        filepattern (str): The pattern to match the data file names.
        temperature (float): The temperature in Kelvin.
        name (str): The name of the TIData.
        detect_equilibrium (bool, optional): Whether to detect equilibrium in the data.
            Defaults to True.
        perWindow (pd.DataFrame, optional): The per-window data.
            Defaults to an empty DataFrame.
        cumulative (pd.DataFrame, optional): The cumulative data.
            Defaults to an empty DataFrame.
        dG (float, optional): The calculated dG value.
            Defaults to 0.
        error (float, optional): The estimated error in the dG value.
            Defaults to 0.
        RT (float, optional): The product of the gas constant and temperature.
            Defaults to 0.
        eqtime (int, optional): The equilibrium time.
            Defaults to 1000.
        upper_walls (float, optional): The upper walls.
            Defaults to None.
        n_lambdas (int, optional): The number of lambdas.
            Defaults to 41.
        lambda_sched (list, optional): The L schedule.
            Defaults to None.
        harmonic_wall (dict, optional): The harmonic wall.
            Defaults to None.
        force_constant (float, optional): The force constant.
            Defaults to 0.
        target_force_constant (float, optional): The target force constant.
            Defaults to 200.
        force_exponent (float, optional): The force exponent.
            Defaults to 6.
        num_steps (int, optional): The number of steps.
            Defaults to 300000.
        data (pd.DataFrame, optional): The data.
            Defaults to an empty DataFrame.

    Methods:
        read(): Reads and processes the colvars data.
        process(): Processes the TIData and calculates dG and error values.
    """

    # pylint: disable=too-many-instance-attributes
    # we might deepen the class heirarchy later, but it works for now.

    eqtime: int = 1000
    upper_walls: float = 0.0
    n_lambdas: int = 41
    lambda_sched: np.ndarray = field(default_factory=np.asarray)
    harmonic_wall: dict = None
    force_constant: int = 0
    target_force_constant: int = 200
    force_exponent: int = 6
    num_steps: int = 300000
    data: pd.DataFrame = field(default_factory=pd.DataFrame)

    def read(self):
        """
        Reads and processes the colvars data from a specified file pattern within
        a given path.

        This method initializes the data processing by reading the first line of
        the file to determine the column names, and then reads the entire file
        into a pandas DataFrame. It filters the data based on the equilibrium time
        (`eqtime`) and adjusts the DataFrame index accordingly.

        Raises:
            FileNotFoundError: If no files match the specified pattern.
            ValueError: If the file is empty or the data format is incorrect.

        Returns:
            None: This method modifies the `data` attribute of the instance in-place
            and does not return anything.
        """
        # Setup and processing of colvars data
        infile = list(self.inpath.glob(self.filepattern))[0]
        with open(infile, "r", encoding="UTF8") as fin:
            first_line = fin.readline()
        columns = re.split(" +", first_line)[1:-1]
        self.data = pd.read_csv(
            infile, delim_whitespace=True, names=columns, comment="#", index_col=0
        )
        self.data = self.data[self.data.index >= self.eqtime][1:]
        self.data.index = self.data.index - self.eqtime

    def process(self):
        """
        Processes the TIData by setting up the lambda schedule if not already set,
        creating a harmonic wall, adjusting the data indices based on the number
        of steps, and finally calculating the thermodynamic integration to determine
        the free energy change (dG) and its error.

        This method first checks if the lambda schedule (`Lsched`) is set, and if
        not, it initializes it with a linear distribution from 1 to 0 over the number
        of lambdas (`n_lambdas`). It then constructs a harmonic wall using the
        specified parameters and adjusts the data indices to correspond to the correct
        lambda values. The method handles edge cases for data indexing and rounds
        the lambda values for better precision. Finally, it processes the thermodynamic
        integration using the adjusted data and updates the per-window and cumulative
        data attributes, along with the calculated dG and error.

        Side effects:
            Modifies several attributes of the instance in-place, including `Lsched`,
            `harmonic_wall`, `data`, `perWindow`, `cumulative`, `dG`, and `error`.

        Returns:
            None
        """
        if not self.lambda_sched:
            self.lambda_sched = np.linspace(1, 0, self.n_lambdas)

        self.harmonic_wall = safep.make_harmonicWall(
            FC=self.force_constant,
            targetFC=self.target_force_constant,
            targetFE=self.force_exponent,
            upperWalls=self.upper_walls,
            targetEQ=self.eqtime,
            numSteps=self.num_steps,
            name=self.name,
            schedule=self.lambda_sched,
        )

        lambdas = (self.data.index.values - 1) // self.harmonic_wall["numSteps"]
        lambdas[0] = 0

        # This is a small hack in case there are extra samples for the last window
        lambdas[lambdas == self.n_lambdas] = self.n_lambdas - 1

        data_lambdas = np.round([self.harmonic_wall["schedule"][i] for i in lambdas], 3)
        self.data.loc[:, "L"] = data_lambdas
        self.data = self.data.iloc[1:]

        self.per_window, self.cumulative = safep.process_TI(
            self.data, self.harmonic_wall, self.lambda_sched
        )
        self.delta_g = np.round(self.cumulative["dG"][1], 1)
        self.error = np.round(self.cumulative["error"][1], 1)


@dataclass
class FEPData(DeltaGData):
    """
    A data class representing FEP (Free Energy Perturbation) data for free energy calculations.

    Attributes:
        inpath (Path): The path to the data files.
        filepattern (str): The pattern to match the data file names.
        temperature (float): The temperature in Kelvin.
        name (str): The name of the FEP data.
        detect_equilibrium (bool, optional): Whether to detect equilibrium in the data.
            Defaults to True.
        perWindow (pd.DataFrame, optional): The per-window data.
            Defaults to an empty DataFrame.
        cumulative (pd.DataFrame, optional): The cumulative data.
            Defaults to an empty DataFrame.
        dG (float, optional): The calculated dG value. Defaults to 0.
        error (float, optional): The estimated error in the dG value.
            Defaults to 0.
        RT (float, optional): The product of the gas constant and temperature.
            Defaults to 0.
        forward (pd.Series, optional): The forward FEP data. Defaults to an empty Series.
        forward_error (pd.Series, optional): The error in the forward FEP data.
            Defaults to an empty Series.
        backward (pd.Series, optional): The backward FEP data.
            Defaults to an empty Series.
        backward_error (pd.Series, optional): The error in the backward FEP data.
            Defaults to an empty Series.

    Methods:
        process(): Process the FEP data and calculate the dG value and error.
        general_plot(width, height, cumulative_ylim, perwindow_ylim):
            Generate a general plot of the FEP data.
        convergence_plot(width, height, fontsize):
            Generate a convergence plot of the FEP data.

    """

    # pylint: disable=too-many-instance-attributes
    # we might deepen the class heirarchy later, but it works for now.

    forward: pd.Series = field(default_factory=pd.Series)
    forward_error: pd.Series = field(default_factory=pd.Series)
    backward: pd.Series = field(default_factory=pd.Series)
    backward_error: pd.Series = field(default_factory=pd.Series)

    def process(self):
        """
        Processes the FEP data to calculate the free energy change (dG) and its error.

        This method computes the gas constant in kcal/(mol K), calculates RT
        (product of the gas constant and temperature), reads and processes FEP
        output files, estimates the free energy using the BAR estimator, and
        calculates convergence data for forward and backward FEP. It updates the
        instance attributes with per-window data, cumulative data, forward and
        backward FEP data, and their respective errors.

        Side effects:
            - Updates several attributes of the instance including `RT`, `perWindow`,
            `cumulative`, `forward`, `forward_error`, `backward`, `backward_error`,
            `dG`, and `error`.

        Returns:
            None
        """
        gas_constant = sp.constants.R / (
            1000 * sp.constants.calorie
        )  # gas constant in kcal/(mol K)
        self.RT = gas_constant * self.temperature  # RT in kcal/mol

        fepout_files = self.inpath.glob(
            self.filepattern
        )  # Resolve any naming regular expressions
        u_nk_site = safep.read_and_process(
            fepout_files,
            self.temperature,
            decorrelate=False,
            detectEQ=self.detect_equilibrium,
        )  # u_nk stores the fep data
        self.per_window, self.cumulative = safep.do_estimation(
            u_nk_site
        )  # Run the BAR estimator on the fep data
        self.forward, self.forward_error, self.backward, self.backward_error = (
            safep.do_convergence(u_nk_site)
        )  # Used later in the convergence plot

        self.delta_g = np.round(self.cumulative.BAR.f.iloc[-1] * self.RT, 1)
        self.error = np.round(self.cumulative.BAR.errors.iloc[-1] * self.RT, 1)

    def general_plot(self, width, height, cumulative_ylim, perwindow_ylim):
        """
        Generates a general plot of the FEP data, including both cumulative and per-window data.

        This method utilizes the `safep.plot_general` function to create a plot
        that visualizes the cumulative and per-window free energy perturbation data.
        The plot is customized with specified dimensions and y-axis limits.

        Args:
            width (int): The width of the plot.
            height (int): The height of the plot.
            cumulative_ylim (tuple): A tuple specifying the y-axis limits for the
                cumulative data plot.
            perwindow_ylim (tuple): A tuple specifying the y-axis limits for the
                per-window data plot.

        Returns:
            tuple: A tuple containing the figure and axes objects from the matplotlib plot.
        """
        fig, axes = safep.plot_general(
            self.cumulative,
            cumulative_ylim,
            self.per_window,
            perwindow_ylim,
            self.RT,
            width,
            height,
            "KDE",
            fontsize=20,
        )

        return fig, axes

    def convergence_plot(self, width, height, fontsize=20):
        """
        Generates a convergence plot for the FEP data, showing forward and backward
        free energy perturbations.

        This method multiplies the forward and backward FEP data and their respective
        errors by the RT value before plotting. It adjusts the plot dimensions
        according to the specified width and height.

        Args:
            width (int): The width of the plot in inches.
            height (int): The height of the plot in inches.
            fontsize (int): Font size used for the plot text. Defaults to 20.

        Returns:
            tuple: A tuple containing the matplotlib figure and axes objects.
        """
        fig, conv_ax = plt.subplots(1, 1)
        conv_ax = safep.convergence_plot(
            conv_ax,
            self.forward * self.RT,
            self.forward_error * self.RT,
            self.backward * self.RT,
            self.backward_error * self.RT,
            fontsize=fontsize,
        )

        fig.set_figwidth(width)
        fig.set_figheight(height)

        return fig, conv_ax
