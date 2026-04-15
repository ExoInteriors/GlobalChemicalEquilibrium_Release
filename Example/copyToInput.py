import os
import shutil
from tools.constants import repo_root


def _case_local_gibbs_wrapper() -> str:
    """Return a small case-local wrapper that runs the repo-level Gibbs.py."""
    return f"""import runpy
import sys
from pathlib import Path

REPO_ROOT = Path(r"{repo_root}")
sys.path.insert(0, str(REPO_ROOT))
runpy.run_path(str(REPO_ROOT / "Gibbs.py"), run_name="__main__")
"""


def copy_inputs(input_dir, version='Sulfur_Version'):
    """
    Copy solver artifacts (solver, param.dat, Gibbs.dat, Gibbs.py) into each generated case folder.
    """
    subfolders = [f.path for f in os.scandir(input_dir) if f.is_dir()]

    # Explicit version sources
    solver_src = os.path.join(repo_root, version, 'solver')
    param_src = os.path.join(repo_root, version, 'param.dat')
    Gibbs_src = os.path.join(repo_root, "Gibbs.dat")
    if not os.path.isfile(Gibbs_src):
        raise FileNotFoundError(f"Missing Gibbs.dat at '{Gibbs_src}'.")
    Gibbs_py_src = os.path.join(repo_root, "Gibbs.py")
    if not os.path.isfile(Gibbs_py_src):
        raise FileNotFoundError(f"Missing Gibbs.py at '{Gibbs_py_src}'.")

    with open(param_src, encoding="utf-8") as pf:
        param_lines = pf.readlines()

    for s in sorted(subfolders):
        print(s)
        shutil.copy(solver_src, os.path.join(s, 'solver'))
        target_param = os.path.join(s, 'param.dat')
        with open(target_param, 'w', encoding="utf-8") as out:
            for line in param_lines:
                # Ensure Gibbs.py path resolves from within each case folder.
                if line.strip().startswith("Gibbs energy file ="):
                    out.write("Gibbs energy file = Gibbs.py\n")
                else:
                    out.write(line)
        shutil.copyfile(Gibbs_src, os.path.join(s, 'Gibbs.dat'))
        with open(os.path.join(s, 'Gibbs.py'), 'w', encoding="utf-8") as out:
            out.write(_case_local_gibbs_wrapper())


if __name__ == "__main__":
    raise SystemExit(
        "Run this via Example/run_gce.py or call copy_inputs(input_dir=..., version=...)."
    )
