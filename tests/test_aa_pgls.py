"""H1 on the dated tree (aa_pgls): the regression against Felsenstein's contrasts, and the species filter.

Five mammals with rounded TimeTree ages (million years): human and chimpanzee split 6.4 My ago, they and the
rhesus macaque 29, mouse and rat 13, primates and rodents 90. Lifespans (AnAge, years) and cysteine content are
made-up but plausible. At lambda = 1 a Brownian-motion GLS slope equals the slope through the origin of the
standardized independent contrasts (Felsenstein 1985), which the test computes by hand on this tree.
"""
import io

import numpy as np
import pandas as pd
from Bio import Phylo

from mitoage_pipelines.longevity.aa_pgls import default_dataset, fit
from mitoage_pipelines.longevity.pgls import shared_path_covariance

NEWICK = "(((human:6.4,chimp:6.4):22.6,macaque:29):61,(mouse:13,rat:13):77);"
ORDER = ["human", "chimp", "macaque", "mouse", "rat"]
LL = dict(human=np.log10(122.5), chimp=np.log10(59.4), macaque=np.log10(40.0), mouse=np.log10(4.0), rat=np.log10(3.8))
CYS = dict(human=0.70, chimp=0.80, macaque=0.95, mouse=1.10, rat=1.05)


def covariance():
    tree = Phylo.read(io.StringIO(NEWICK), "newick")
    tips = {t.name: i for i, t in enumerate(tree.root.get_terminals())}
    return shared_path_covariance(tree.root, tips)[np.ix_([tips[s] for s in ORDER], [tips[s] for s in ORDER])]


def contrasts(v):
    """Felsenstein's standardized contrasts on NEWICK, by hand: (x1 - x2) / sqrt(v1 + v2) at each node; the
    node's value is the inverse-variance mean and its branch grows by v1 v2 / (v1 + v2)."""
    def node(a, va, b, vb):
        return (a - b) / np.sqrt(va + vb), (a / va + b / vb) / (1 / va + 1 / vb), va * vb / (va + vb)
    c1, hc, vhc = node(v["human"], 6.4, v["chimp"], 6.4)
    c2, p, vp = node(hc, 22.6 + vhc, v["macaque"], 29.0)
    c3, r, vr = node(v["mouse"], 13.0, v["rat"], 13.0)
    c4, _, _ = node(p, 61.0 + vp, r, 77.0 + vr)
    return np.array([c1, c2, c3, c4])


def test_gls_at_lambda_one_is_the_contrast_regression():
    x, y = contrasts(LL), contrasts(CYS)
    expected = (x @ y) / (x @ x)                       # through the origin
    got = fit(covariance(), np.array([CYS[s] for s in ORDER]), np.array([LL[s] for s in ORDER]), lam=1.0)
    assert np.isclose(got["slope"], expected, rtol=1e-10)
    assert np.isclose(expected, -0.28735, atol=5e-6)  # the number itself, so that a change of the case is seen
    r = y - expected * x                              # the slope's standard error from the contrasts (df = 4 - 1)
    assert np.isclose(got["se"], np.sqrt((r @ r / 3) / (x @ x)), rtol=1e-10)


def test_lambda_zero_is_ordinary_least_squares():
    C = covariance(); y = np.array([CYS[s] for s in ORDER]); x = np.array([LL[s] for s in ORDER])
    assert np.allclose(np.diag(C), 90.0)              # ultrametric: at lambda 0 every species weighs the same
    assert np.isclose(fit(C, y, x, lam=0.0)["slope"], np.polyfit(x, y, 1)[0], rtol=1e-10)


def test_default_dataset():
    sp = pd.DataFrame({
        "record_status": ["original", "new", "retired", "updated", "renamed", "new", "new"],
        "lifespan_quality": ["high", "acceptable", "high", "questionable", "low", "acceptable", None],
        "lifespan_record_type": ["captive", "wild", "captive", "observed", "estimated", "estimated", "wild"],
    })
    assert default_dataset(sp).tolist() == [True, True, False, False, False, False, False]
