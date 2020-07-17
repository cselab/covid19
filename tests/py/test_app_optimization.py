import shutil
import subprocess
import tempfile
import unittest
import os

class TestAppEvidence(unittest.TestCase):
    def setUp(self):
        """Create a temporary directory before each test."""
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Delete the temporary directory after each test."""
        shutil.rmtree(self.tmp_dir)

    def test_optimize_single(self):
        """Test that optimize.py runs without crashing."""
        cmd = [
            os.path.normpath(os.path.join(__file__, '..', '..', '..', 'applications', 'optimization', 'optimize.py')),
            '--compModel', 'country.reparam.sird_dint_nogamma.nbin',
            '--country', 'switzerland',
            '--silent',
            '--silentPlot',
            '--nSamples', '8',
            '--nGenerations', '2',
            '--useInfections',
            '--useDeaths',
            '--dataFolder', self.tmp_dir,
        ]

        subprocess.check_call(cmd, cwd=self.tmp_dir, stdout=subprocess.DEVNULL)

