import re
from pathlib import Path

import numpy as np


def test_readme_quick_start_example_runs():
    readme = Path(__file__).parents[1] / "README.md"
    content = readme.read_text(encoding="utf-8")
    quick_start = content.split("## Quick Start\n", maxsplit=1)[1].split(
        "\n## ", maxsplit=1
    )[0]
    code_blocks = re.findall(r"```python\s*\n(.*?)```", quick_start, re.DOTALL)

    assert len(code_blocks) == 1

    namespace = {"__name__": "__main__"}
    exec(compile(code_blocks[0], "README.md", "exec"), namespace)

    results = namespace["res"]
    assert np.isfinite(results["CplsP"].data).all()
    assert np.isfinite(results["NplsP"].data).all()
