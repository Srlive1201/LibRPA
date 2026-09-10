"""Stage an existing Si dataset without modifying its numerical files."""
import argparse
from pathlib import Path
import re
import shutil

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("dataset", type=Path)
parser.add_argument("destination", type=Path)
parser.add_argument("--legacy", action="store_true", help="Run the original PR revision")
args = parser.parse_args()
source = args.dataset.resolve(strict=True)
args.destination.mkdir(parents=True, exist_ok=False)
target = args.destination / "dataset"
target.mkdir()
manifests = {"qsgw_input.contract", "vxc_scf_manifest", "vxc_band_manifest"}
for item in source.iterdir():
    if item.name in manifests:
        text = item.read_text()
        if not args.legacy:
            text = text.replace("librpa-qsgw-input-contract-v1", "librpa-qsgw-input-contract-v2")
            text = text.replace("librpa-qsgw-vxc-manifest-v2", "librpa-qsgw-vxc-manifest-v3")
            text = re.sub(r"(?<!\S)[0-9a-fA-F]{64}\s+", "", text)
            text = re.sub(r"\bsha256\s+", "", text)
        (target / item.name).write_text(text)
    else:
        (target / item.name).symlink_to(item, target_is_directory=item.is_dir())
run = args.destination / "librpa"
run.mkdir()
shutil.copyfile(Path(__file__).with_name("librpa.in"), run / "librpa.in")
print(run.resolve())
