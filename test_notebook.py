from testbook import testbook

OLD_NOTEBOOK_PATH =  'SAFEP_Tutorial_Notebook.ipynb'
NEW_NOTEBOOK_PATH = 'SAFEP_Tutorial_Notebook_New.ipynb'

def test_delta_g_values():
    with testbook(OLD_NOTEBOOK_PATH, execute=range(0, 11)) as tb, testbook(NEW_NOTEBOOK_PATH, execute=True) as tc:
        tp1 = tb.ref("site")
        tp2 = tc.ref("delta_g")
        assert tp1.delta_g.item() == tp2
