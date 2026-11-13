"""Import script modules without introducing a Python package/restructuring the CLI."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
for folder in (ROOT / 'scripts' / 'xplora', ROOT / 'scripts' / 'iraman'):
    sys.path.insert(0, str(folder))
