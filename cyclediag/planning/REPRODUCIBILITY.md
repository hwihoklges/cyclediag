# Bounded reproducibility metadata

Implemented entry points:

- `cyclediag.api.extract_features`: JSON-safe `DataFrame.attrs['provenance']`,
  plus normalized unit metadata/warnings. Caller configuration is deep-copied;
  cycle iterables are materialized once.
- `cyclediag.api.diagnose_csv`: explicit result `provenance`, including extraction
  configuration digest and diagnosis options. Feature attrs retain extraction provenance.
- CLI `extract`: all branches (basic, both LGES selectors, CV-only) attach metadata;
  CSV output uses the checked sidecar helper. Parquet retains its existing writer;
  no portable sidecar/integrity guarantee is made for Parquet.
- CLI `diagnose`: only its feature CSV has a checked sidecar, recording the
  end-to-end result provenance. Scores/screens remain ordinary CSVs.

## Content and privacy

Provenance records an explicit algorithm schema label, SHA-256 input/config/code
digests, available Git commit/dirty status, and a small allowlist of installed
package versions (CycleDiag, numpy, pandas, scipy, scikit-learn, pyarrow).
No timestamps, hostname, environment variables, absolute paths, raw configuration
values or cell identifiers are added to provenance. Configuration hashes include
dataclass defaults, ColumnMap/units, cycle selections and entry-point selectors;
resolved dQ/dV defaults and external diagnosis-config bytes are included too.
Derived input filenames/cell IDs can affect the configuration digest without being
disclosed. These hashes are not anonymization guarantees for low-entropy identifiers.

File inputs hash original bytes in streaming blocks. DataFrame inputs instead hash
the normalized table (ordered values/index, column/dtype schema and unit attrs),
using chunked pandas row hashes under SHA-256. These source kinds are explicitly
different and are not interchangeable; DataFrame hashes are pandas-version-bound,
not a cryptographic canonical encoding of every possible Python object. Inputs
must remain unchanged during extraction; concurrent modification is not locked.

Code SHA-256 covers package Python/JSON/TOML files (relative paths and content),
excluding tests/caches. It captures local scientific/configuration edits even for
dirty trees, but is not a complete build artifact or environment lock. Git uses
argument arrays, no shell and a three-second timeout per command. A commit is only
attributed when the exact package API and initializer are tracked under that root;
an unrelated parent checkout around an installed package is not sufficient.
Unavailable Git returns null commit/dirty values. Dirty status is repository-wide.

## Portable feature CSVs

Use `cyclediag.provenance.save_features_csv(frame, path, source_path=...)` and
`read_feature_csv(path)`. The fixed companion name is the CSV filename followed by
`.metadata.json`; no paths supplied by a manifest are followed. Loading requires
the sidecar and rejects changed CSV bytes or incompatible schema. This detects
accidental corruption/stale pairs, not malicious replacement of both files.

The sidecar preserves allowlisted attrs: provenance, unit metadata/schema/warnings,
missing unit fields, feature units and scientific warnings. Other arbitrary attrs
are intentionally omitted. List/dict/tuple cells (including warning and feature-unit
columns) use JSON and are reconstructed without eval/pickle; nonfinite JSON values
become null, tuples become lists. CSV scalar dtype inference remains ordinary pandas
behavior; this is not a lossless DataFrame/index serializer. Existing feature cell
IDs, paths, warnings and source column labels are **not redacted** by this export.
Review exported data/metadata before sharing it.

Both destinations reject known source collisions (explicit `source_path`, or the
feature `file` column), symlinks and existing hard-link aliases. When a table has no
source information, callers must pass `source_path` to protect that input. Existing
output files may be replaced. Temporary files are independently created in the
destination directory and cleaned up; each replacement is atomic, but the two-file
pair is not transactional. A crash between replacements produces a detectable
integrity failure, not silently trusted metadata.

## Remaining scope

Low-level extraction/diagnosis functions, `diagnose_dataframe`, folder/batch/DOE,
peak pipelines, training/prediction and standalone exporters are unchanged. They
do not automatically construct provenance. Pandas operations can discard attrs;
retain the explicit `diagnose_csv` result or use the checked feature CSV helpers.
Full dependency locking, execution/environment capture, signed manifests, complete
pipeline provenance and controlled real-data scientific validation remain roadmap
items. Metadata does not establish scientific validity or exact cross-platform
numerical reproducibility. No scientific algorithms are changed by this support.