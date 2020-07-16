import subprocess
import unittest


class TestCountrySIR(unittest.TestCase):
    def test_sir_vs_seir(self):
        """Test comparison of SIR and SEIR with interventions."""
        cwd = "../../applications/gui_paper/doc"
        name = "solve_sir_vs_seir.py"
        subprocess.check_call("./{} --test".format(name), shell=True, cwd=cwd)
