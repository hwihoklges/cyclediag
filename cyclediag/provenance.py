"""Bounded provenance and checked, portable feature CSVs (no environment capture)."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import hashlib
from importlib.metadata import PackageNotFoundError, version
import io
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd

ALGORITHM_SCHEMA = "cyclediag_scientific_v2"
ATTR_KEYS = ("unit_metadata", "unit_warnings", "unit_schema", "unit_missing_fields",
             "scientific_warnings", "feature_units", "provenance")


def json_safe(value):
    """Strict value serialization: never introspect arbitrary objects or repr them."""
    if is_dataclass(value) and not isinstance(value, type):
        return json_safe(asdict(value))
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        if not all(isinstance(k, str) for k in value):
            raise TypeError("Metadata keys must be strings")
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    raise TypeError("Unsupported metadata value type")


def canonical_json(value) -> str:
    return json.dumps(json_safe(value), sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)


def config_hash(config) -> str:
    return hashlib.sha256(canonical_json(config).encode("utf-8")).hexdigest()


def file_hash(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def frame_metadata(frame):
    return json_safe({k: frame.attrs[k] for k in ATTR_KEYS if k in frame.attrs})


def dataframe_hash(frame) -> str:
    """Hash normalized values/index in row order, schema and unit attrs; pandas-version bound."""
    schema = {"columns": list(frame.columns), "dtypes": [str(t) for t in frame.dtypes],
              "index_dtype": str(frame.index.dtype), "index_names": list(frame.index.names),
              "units": {k: frame.attrs[k] for k in ATTR_KEYS[:4] if k in frame.attrs}}
    digest = hashlib.sha256(canonical_json(schema).encode("utf-8"))
    for start in range(0, len(frame), 10000):
        hashed = pd.util.hash_pandas_object(frame.iloc[start:start + 10000], index=True)
        digest.update(hashed.to_numpy(dtype="<u8").tobytes())
    return digest.hexdigest()


def code_metadata(package_dir=None):
    package = Path(package_dir or Path(__file__).resolve().parent).resolve()
    digest = hashlib.sha256()
    for path in sorted(package.rglob("*")):
        relative = path.relative_to(package)
        if (path.is_file() and path.suffix in (".py", ".json", ".toml")
                and not any(p in ("tests", "__pycache__", ".pytest_cache") for p in relative.parts)):
            digest.update(relative.as_posix().encode("utf-8") + b"\0")
            digest.update(file_hash(path).encode("ascii"))
    result = {"source_sha256": digest.hexdigest(), "commit": None, "dirty": None}

    def git(*args):
        return subprocess.run(["git", "-C", str(package), *args], check=True,
                              capture_output=True, text=True, timeout=3, shell=False).stdout.strip()

    try:
        root = Path(git("rev-parse", "--show-toplevel")).resolve()
        relative = package.relative_to(root)
        # A parent repository is not evidence: this exact installed code must be tracked.
        tracked = git("ls-files", "--full-name", "--", "__init__.py", "api.py").splitlines()
        required = {(relative / name).as_posix() for name in ("__init__.py", "api.py")}
        if required.issubset(set(tracked)):
            result.update(commit=git("rev-parse", "HEAD"),
                          dirty=bool(git("status", "--porcelain", "--untracked-files=normal")))
    except (OSError, subprocess.SubprocessError, ValueError):
        pass
    return result


def build_provenance(source, normalized, *, config, column_map, cycles=None, selectors=None):
    effective = {"config": json_safe(config), "column_map": json_safe(column_map),
                 "cycles": cycles, "selectors": selectors or {},
                 "unit_schema": normalized.attrs.get("unit_schema"),
                 "unit_metadata": normalized.attrs.get("unit_metadata", {})}
    if hasattr(config, "resolved_dqdv_config"):
        effective["resolved_dqdv_config"] = json_safe(config.resolved_dqdv_config())
    external = getattr(config, "diagnosis_config_path", None)
    if external:
        effective["diagnosis_config_sha256"] = file_hash(external)
    packages = {}
    for name in ("cyclediag", "numpy", "pandas", "scipy", "scikit-learn", "pyarrow"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    is_frame = isinstance(source, pd.DataFrame)
    return {"schema_version": 1, "algorithm_schema": ALGORITHM_SCHEMA,
            "input": {"kind": "normalized_dataframe_pandas_v1" if is_frame else "file_bytes",
                      "sha256": dataframe_hash(normalized) if is_frame else file_hash(source)},
            "config_sha256": config_hash(effective), "code": code_metadata(),
            "packages": packages}


def attach_provenance(table, normalized, provenance):
    table.attrs.update(frame_metadata(normalized))
    table.attrs["provenance"] = json_safe(provenance)
    return table


def save_features_csv(frame, path, *, source_path=None):
    """Save CSV + fixed .metadata.json sidecar. Reject known input paths and symlinks.

    Only allowlisted attrs are exported. Structured cells use JSON, not Python repr.
    Each replacement is atomic, not the pair; an interrupted pair fails verification.
    """
    target = Path(path).absolute()
    sidecar = target.with_name(target.name + ".metadata.json")
    sources = [source_path] if source_path is not None else []
    if "file" in frame:
        sources.extend(v for v in frame["file"].dropna().unique() if isinstance(v, str) and v)
    for dest in (target, sidecar):
        if dest.is_symlink() or any(
            dest.resolve() == Path(src).resolve() or
            (dest.exists() and Path(src).exists() and os.path.samefile(dest, src))
            for src in sources
        ):
            raise ValueError("Feature export must not overwrite a source or symlink")
    if not frame.columns.is_unique or not all(isinstance(c, str) for c in frame.columns):
        raise ValueError("Feature CSV requires unique string columns")
    data = frame.copy(deep=True)
    json_columns = []
    for col in data.columns:
        if any(isinstance(v, (dict, list, tuple)) for v in data[col]):
            json_columns.append(col)
            data[col] = data[col].map(canonical_json)
    attrs = frame_metadata(frame)
    target.parent.mkdir(parents=True, exist_ok=True)
    temps = []
    try:
        for _ in range(2):
            fd, name = tempfile.mkstemp(dir=target.parent, prefix=".cyclediag-", suffix=".tmp")
            os.close(fd)
            temps.append(Path(name))
        # Object conversion avoids numpy's NaN-to-string RuntimeWarning on older pandas.
        data = data.astype(object).where(pd.notna(data), "")
        data.to_csv(temps[0], index=False, encoding="utf-8")
        manifest = {"schema_version": 1, "csv_sha256": file_hash(temps[0]),
                    "columns": list(frame.columns), "json_columns": json_columns, "attrs": attrs}
        temps[1].write_text(canonical_json(manifest), encoding="utf-8")
        os.replace(temps[0], target)
        os.replace(temps[1], sidecar)
    finally:
        for temp in temps:
            temp.unlink(missing_ok=True)
    return sidecar


def read_feature_csv(path):
    """Require matching sidecar; never evaluate Python cell strings or follow manifest paths."""
    path = Path(path)
    manifest = json.loads(path.with_name(path.name + ".metadata.json").read_text(encoding="utf-8"))
    payload = path.read_bytes()
    if manifest.get("schema_version") != 1 or hashlib.sha256(payload).hexdigest() != manifest.get("csv_sha256"):
        raise ValueError("Feature CSV integrity/schema mismatch")
    frame = pd.read_csv(io.BytesIO(payload), converters={
        col: json.loads for col in manifest["json_columns"]
    }) if manifest["columns"] else pd.DataFrame()
    if list(frame.columns) != manifest["columns"]:
        raise ValueError("Feature CSV column mismatch")
    frame.attrs.update(json_safe(manifest["attrs"]))
    return frame