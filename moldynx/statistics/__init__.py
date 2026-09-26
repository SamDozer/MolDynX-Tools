"""Statistics: descriptive stats, time-series helpers, correlation, bootstrap CIs."""

from moldynx.statistics.descriptive import describe, summary_frame, bootstrap_ci  # noqa: F401
from moldynx.statistics.timeseries import (  # noqa: F401
    moving_average, rolling_std, running_mean, block_average, block_average_sem,
    plateau_detection,
)
from moldynx.statistics.correlation import correlation_matrices  # noqa: F401
from moldynx.statistics.autocorr import (  # noqa: F401
    statistical_inefficiency, describe_correlated, detect_equilibration, drift,
)
