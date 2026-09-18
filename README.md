# CycleDiag

Standalone **cycle / voltage-profile diagnosis** package.

GUI 없음 · `pne_studio2` 불필요 · 단독 구동.

**Version:** 1.0.0

## 2026-09 정확성 리뷰

[과학적 계약·검증 계획·단계별 로드맵](cyclediag/planning/SCIENTIFIC_REVIEW_2026-09.md)을 먼저 확인하세요. 내부 단위는 A/Ah/V/s, `f_Q_spec`은 mAh/g입니다. 명시적 `ColumnMap.units`와 헤더 단위의 충돌은 오류이며, 단위 없는 입력은 canonical 가정 경고를 남깁니다. `DataFrame.attrs`만으로 영속 추적성을 보장하지 않습니다.

**unknown ≠ 건강함.** 미교정 점수는 확률이 아니며 full-cell LLI/LAM/kinetics proxy는 원인을 유일하게 식별하지 못합니다. 구현된 수치 방어와 외부 과학 검증 완료를 구분합니다. 기존 연구 로드맵은 보존하고 새 문서에 단계별 의존성·인수 기준을 연결했습니다.

## Example fixtures (Git LFS) — **DOE1 / DOE2 / DOE3**

비교할 때 **`DOE1`**, **`DOE2`**, **`DOE3`** 만 말하면 됩니다.  
상세: [`example/fixtures/README.md`](example/fixtures/README.md) · [`manifest.json`](example/fixtures/manifest.json)

| ID | 비교 | Arms |
|----|------|------|
| **DOE1** | SJ900 wet vs dry | set1 · set4 |
| **DOE2** | SJ900 dry vs SJ1300 dry | SJ1300_dry (+ DOE1 set4 ref) |
| **DOE3** | **양극** Bimodal vs S83S | Ch109–111 · Ch103–105 |
| halfcell | BOL OCP C/20 | anode · cathode (**aged 없음**) |

```bash
git lfs install
git lfs pull
python run_cyclediag.py extract --input example/fixtures/doe/DOE1/set4_SJ900/M01Ch025_raw.csv --out /tmp/f.csv
```
## Cursor Cloud Agents

이 레포를 Cloud Agents에 연결하면 `.cursor/environment.json`의 install이 의존성을 설치합니다.

1. [Cloud Agents dashboard](https://cursor.com/dashboard/cloud-agents#environments) → GitHub에서 이 레포 선택
2. Environment 셋업(에이전트 자동 설치 권장)
3. 태스크 예: `DOE1`의 `example/fixtures/doe/DOE1/set4_SJ900/M01Ch025_raw.csv`로 diagnose 돌리고 결과 요약해 PR 열어줘

## Install

```bash
pip install -r requirements.txt
pip install -e "./cyclediag[dev]"
```

## Library

```python
from cyclediag import diagnose_csv, extract_features
from cyclediag.models.predict import predict_features

result = diagnose_csv("cell_raw.csv")
print(result["scored"][["cycle", "SoHQ", "anomaly_score", "flag"]].head())
```

## CLI

```bash
python -m cyclediag --help
python run_cyclediag.py extract --input raw.csv --out features.csv
python run_cyclediag.py diagnose --input raw.csv --out-dir out/diag
python run_cyclediag.py predict --features features.csv --out scores.csv
python run_cyclediag.py report --input-dir path/to/folder
python run_cyclediag.py compare-doe --doe DOE2 --out example/output/DOE2_compare
python run_cyclediag.py peaks export --input raw.csv --out-dir example/docs/features --cell-id Cell01
```

## SoHQ BP + voltage-profile presentation (tagged)

Tagged-routine SoHQ regimes (BP1/BP2, 3 zones) and every-N-cycle V–Q / dQ/dV around each BP:

```bash
# all default cells (set4 SJ900 + SJ1300 dry)
python -m cyclediag.tools.run_sohq_bp_presentation --out example/output/sohq_bp_presentation

# smoke (one cell, coarser stride)
python -m cyclediag.tools.run_sohq_bp_presentation --cells M01Ch022 --step 40 --out /tmp/bp_pres
```

Outputs per cell: `*_sohq_dsohq_regimes.png`, `*_BP1_BP2_VQ_dQdV.png`, regime slope CSV.

**Cloud / CI:** `cyclediag-ci`는 Python 3.10/3.14, Ubuntu/Windows의 3개 조합에서 전체 단위 테스트를 **LFS 없이** 실행합니다. fixture skip은 JUnit과 `-rs`로 공개합니다. Presentation은 별도 수동 `cyclediag-integration` workflow에서 `run_integration=true`를 선택해야 하며 M01Ch022 raw 한 파일만 받습니다. 의존성 목록과 결과를 artifact로 보존합니다.
Cursor Cloud Agents: clone this repo → install via `.cursor/environment.json` → run the same command.

## dVdQ@SOC0 Origin OLE PowerPoint

HPPC와 같이 Origin 그래프를 PowerPoint에 OLE로 붙여 넣습니다. Origin 2025 + 데스크톱 PowerPoint가 필요합니다.

```bash
pip install -e "./cyclediag[origin]"
python -m cyclediag.origin_ppt
python -m cyclediag.origin_ppt --skip-profiles
python cyclediag/tools/build_dvdq_soc0_slides.py --skip-raw
```

출력: `example/output/dvdq_soc0_slides/`

## Bounded reproducibility

High-level `extract_features` / `diagnose_csv` and CLI extraction now record input,
configuration and code digests. Feature CSV exports support checked JSON sidecars
for provenance, units and scientific warnings. See
[reproducibility support and limits](cyclediag/planning/REPRODUCIBILITY.md)
for exact integrated entry points, privacy behavior and remaining roadmap scope.

## Layout

```
cyclediag/          # package (diagnosis params under diagnosis/config/)
  origin_ppt/       # dVdQ@SOC0 Origin OLE deck
example/fixtures/   # DOE raw.csv + halfcell/ BOL OCP (Git LFS)
.cursor/            # Cloud Agent environment
run_cyclediag.py
run_export_cycle_indicators.py
```

`pne_studio2` Diagnosis 탭 · rest_voltage는 이 패키지를 사용합니다.
