"""rcmm - Python side of rates-curve-nlp.

A numpy port of the header-only C++ curve engine (``include/rcmm/curve.hpp``) and factor model
(``factors.hpp``), used as the independent reference for the C++ known-answer tests; loaders for the
downloaded data and for every ``results/*.json`` the pipeline writes; and the small statistics toolkit
(bootstrap / permutation / block-bootstrap CIs) shared by the text, hedge and event scripts.

    >>> from rcmm import curve, data
    >>> rows = data.read_cmt("data/sample/treasury_cmt_sample.csv")
    >>> inst = data.cmt_instruments(rows.loc["2025-04-14"])
    >>> b = curve.Bootstrap().fit(inst); kr = curve.KernelRidge().fit(inst, 1e-4)
    >>> round(b.par_yield(10.0), 6), round(kr.max_error_bp(inst), 3)
    (0.0438, 0.092)
"""
from . import curve, data, factors, results, stats
from .curve import (Bootstrap, InstKind, Instrument, KernelRidge, MonotoneConvex, StepFrontEnd,
                    forward_roughness, key_rate_dv01, price_instrument)
from .data import CMT_NAMES, CMT_T, cmt_instruments, read_cmt
from .factors import FactorModel

__all__ = ["curve", "data", "factors", "results", "stats", "Bootstrap", "InstKind", "Instrument", "KernelRidge",
           "MonotoneConvex", "StepFrontEnd", "forward_roughness", "key_rate_dv01", "price_instrument", "CMT_NAMES",
           "CMT_T", "cmt_instruments", "read_cmt", "FactorModel"]
__version__ = "0.1.0"
