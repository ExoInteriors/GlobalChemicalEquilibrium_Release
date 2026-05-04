"""Configure Matplotlib for non-interactive workflow runs."""

import os
import tempfile


def configure_matplotlib_cache() -> None:
    if "MPLCONFIGDIR" in os.environ:
        return
    cache_dir = os.path.join(tempfile.gettempdir(), "gce_matplotlib")
    os.makedirs(cache_dir, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = cache_dir
