from pathlib import Path
import runpy


SCRIPT_PATH = Path(__file__).resolve().parent / "scripts" / "prep" / "prepare_GSE115978.py"


if __name__ == "__main__":
    runpy.run_path(str(SCRIPT_PATH), run_name="__main__")
