import shutil
import subprocess
import tempfile
import unittest

class TestAppEvidence(unittest.TestCase):
    def setUp(self):
        """Create a temporary directory before each test."""
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Delete the temporary directory after each test."""
        shutil.rmtree(self.tmp_dir)

    def test_sample_knested(self):
        """Test that sample_knested.py runs without crashing."""
        cmd = [
            './sample_knested.py',
            '--compModel', 'country.reparam.sird_dint_nogamma.nbin',
            '--country', 'switzerland',
            '--silent',
            '--silentPlot',
            '--nSamples', '20',
            '--useInfections',
            '--useDeaths',
            '--dataFolder', self.tmp_dir,
            '--test',
        ]

        subprocess.check_call(cmd, cwd='../../applications/evidence', stdout=subprocess.DEVNULL)
