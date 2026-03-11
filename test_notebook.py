from testbook import testbook
import numpy as np

tb1 = @testbook('SAFEP_Tutorial_Notebook_Old.ipynb', execute=range(0, 11))

    
tb2 = @testbook('SAFEP_Tutorial_Notebook.ipynb', execute=True)


def test_stdout(tb1, tb2):
    site = tb1.ref("site") 
    site2 = tb1.ref("site")
    assert site.delta_g.item() == site2