from unittest import TestCase

class TestCaseEx(TestCase):
    def assertRelative(self, a, b, tolerance):
        relative = abs((a - b) / abs(a))
        if relative > tolerance:
            self.fail(f"assertRelative failed: |{a} - {b}| relative "
                      f"error {relative} larger than tolerance {tolerance}.")
