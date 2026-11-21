"""Read the historical software catalogue and reconstruct a named release.

Original ZIP hashes identify provenance. Member hashes identify extracted bytes.
A single object may occur under several names or in several releases.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

ARCHIVE = Path(__file__).resolve().parents[1] / "archive"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


class Catalogue:
    def __init__(self):
        self.handles = [ZipFile(ARCHIVE / "sources.zip")]
        try:
            self.handles.extend(ZipFile(p) for p in sorted((ARCHIVE / "support").glob("*.zip")))
            self.origins = json.loads(self.handles[0].read("manifest.json"))["origins"]
            self.objects = {}
            for handle in self.handles:
                for name in handle.namelist():
                    if name.startswith("objects/"):
                        digest = name.removeprefix("objects/")
                        if digest in self.objects:
                            raise ValueError(f"Duplicate physical object {digest}")
                        self.objects[digest] = handle
        except BaseException:
            self.close()
            raise

    def close(self):
        for handle in self.handles:
            handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def data(self, record):
        digest = record["sha256"]
        value = self.objects[digest].read("objects/" + digest)
        if sha256(value) != digest or len(value) != record["bytes"]:
            raise ValueError(f"Damaged archive object {digest}")
        return value

    def extract(self, name, directory):
        """Restore source and support files to a new root. Never run archived code here."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=False)
        for member, record in self.origins[name]["members"].items():
            path = PurePosixPath(member)
            if path.is_absolute() or ".." in path.parts or "\\" in member:
                raise ValueError(f"Unsafe member path {member}")
            target = directory / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(self.data(record))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", nargs=2, metavar=("ORIGIN", "NEW_DIRECTORY"))
    args = parser.parse_args()
    with Catalogue() as catalogue:
        if args.extract:
            catalogue.extract(*args.extract)
        else:
            for name, origin in catalogue.origins.items():
                print(f"{name}  ({len(origin['members'])} members)")
