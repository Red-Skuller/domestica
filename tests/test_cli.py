import pytest
from domestica.io_utils import clean_seq

def test_clean_seq():
    """Simple test to ensure amino acid cleaning works"""
    raw_seq = "M-A*C T"
    assert clean_seq(raw_seq) == "MACT"