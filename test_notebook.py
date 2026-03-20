from testbook import testbook
from matplotlib.testing.compare import compare_images
from pathlib import Path

OLD_NOTEBOOK_PATH =  'SAFEP_Tutorial_Notebook.ipynb'
NEW_NOTEBOOK_PATH = 'SAFEP_Tutorial_Notebook_New.ipynb'
ROOT = "."
TEST_FILES = Path(f"{ROOT}/test_output_notebook_plots")
BOUND_FEP_PATH = Path(f"{ROOT}/stepB_alchemy_site/sample_outputs/")
RESTRAINT_PERTURBATION_PATH = Path(f"{ROOT}/stepC_restraint_perturbation/sample_outputs/")
BULK_FEP_PATH = Path(f"{ROOT}/stepD_alchemy_bulk/sample_outputs/")

@pytest.fixture(scope="session")
def shared_data_folder() -> Path:
    # Create a directory named "shared_data" in the base temporary directory
    data_dir: Path = TEST_FILES.mktemp("shared_data")
    yield data_dir
    
def test_all_notebook_values():
    """
    Reads in values from the old and new 
    testing notebook, then compares equivalent 
    output values.

    """
    with testbook(OLD_NOTEBOOK_PATH, execute=range(0, 12)) as tb, testbook(NEW_NOTEBOOK_PATH, execute=True) as tc:
        # test if deltaG is equivalent
        tp1 = tb.ref("site")
        tp2 = tc.ref("bound_delta_g")
        assert tp1.delta_g.item() == tp2

        tp1 = tb.ref("DBC")
        tp2 = tc.ref("dbc_delta_g")
        assert tp1.delta_g.item() == tp2

def test_general_figure_plot_equivalence():
    """
    Test if the bound_generalFigure.pdf figures are equivalent.
    
    """
    bound_gen_fig = "bound_generalFigures.pdf"
    assert None == compare_images(str(BOUND_FEP_PATH.joinpath(bound_gen_fig)), str(TEST_FILES.joinpath(bound_gen_fig)), tol=0.01)
    
def test_convergence_plot_equivalence():
    """
    Test if bound_convergence.pdf figures are equivalent.

    """
    bound_gen_fig = "bound_convergence.pdf"
    assert None == compare_images(str(BOUND_FEP_PATH.joinpath(bound_gen_fig)), str(TEST_FILES.joinpath(bound_gen_fig)), tol=0.01)

def test_dbc_plot_equivalence():
    """
    Test if TI_general.pdf figures are equivalent.

    """
    bound_gen_fig = "TI_general.pdf"
    assert None == compare_images(str(RESTRAINT_PERTURBATION_PATH.joinpath(bound_gen_fig)), str(TEST_FILES.joinpath(bound_gen_fig)), tol=0.01)
