"""
Tests for the Subject Class 
"""
from cohort import Subject
import pytest

class TestSubject(Subject):
    """
    Class to test Subject Class.
    """
    def setUp(self):
        """
        Initialization class
        """
        self.subject = Subject()
    def test_parameters(self):
        """
        Test if all parameters are within expected probability distribution.
        """
