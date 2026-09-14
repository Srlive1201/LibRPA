"""Stage an existing Si dataset without modifying its numerical files."""
import argparse
from pathlib import Path
import shutil

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("dataset", type=Path)
parser.add_argument("destination", type=Path)
args = parser.parse_args()
source = args.dataset.resolve(strict=True)
args.destination.mkdir(parents=True, exist_ok=False)
target = args.destination / "dataset"
target.mkdir()
for item in source.iterdir():
    (target / item.name).symlink_to(item, target_is_directory=item.is_dir())
run = args.destination / "librpa"
run.mkdir()
shutil.copyfile(Path(__file__).with_name("librpa.in"), run / "librpa.in")
print(run.resolve())
