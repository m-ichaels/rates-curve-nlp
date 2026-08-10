import json
import os

import pytest

from rcmm import data as D

ROOT = D.ROOT
SAMPLE = os.path.join(ROOT, "data", "sample", "treasury_cmt_sample.csv")
RESULTS = os.path.join(ROOT, "results")


@pytest.fixture(scope="session")
def cmt():
    return D.read_cmt(SAMPLE)


@pytest.fixture(scope="session")
def inst(cmt):
    return D.cmt_instruments(cmt.loc["2025-04-14"])


@pytest.fixture(scope="session")
def cpp_curve():
    """The C++ `rcmm curve --date 2025-04-14` output checked in under results/."""
    with open(os.path.join(RESULTS, "curve_2025-04-14.json")) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def cpp_factors():
    with open(os.path.join(RESULTS, "factors.json")) as f:
        return json.load(f)
