from unittest import TestCase

class TestCaseEx(TestCase):
    def assertRelative(self, a, b, tolerance):
        relative = abs((a - b) / abs(a))
        if relative > tolerance:
            self.fail(f"assertRelative failed: |{a} - {b}| not under "
                      f"relative tolerance of {tolerance}.")
