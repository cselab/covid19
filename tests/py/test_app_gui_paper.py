import subprocess
import unittest

class TestAppGUIPaper(unittest.TestCase):
    def test_sir_vs_seir(self):
        """Test comparison of SIR and SEIR with interventions."""
        cwd = '../../applications/gui_paper/doc'
        subprocess.check_call(['./solve_sir_vs_seir.py', '--test'], cwd=cwd)
