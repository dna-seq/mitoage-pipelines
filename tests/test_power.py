"""H2's power: the detectable effect is se times (two-sided critical value + the power's quantile)."""
import numpy as np

from mitoage_pipelines.longevity.power import detectable


def test_textbook_case():
    # large df: the normal quantiles 1.95996 (alpha 0.05, two-sided) and 0.84162 (power 0.8)
    assert np.isclose(detectable(1.0, 0.05, 10**7), 1.95996 + 0.84162, atol=1e-4)
    assert np.isclose(detectable(0.03, 0.05, 10**7), 0.03 * 2.80158, atol=1e-5)


def test_bh_threshold_of_one_discovery_in_mammals():
    # alpha = 0.05 / 5259: two-sided normal critical value 4.42808 (qnorm(1 - 0.05 / 5259 / 2))
    assert np.isclose(detectable(1.0, 0.05 / 5259, 10**7), 4.42808 + 0.84162, atol=1e-4)
