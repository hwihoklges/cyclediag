# dQ/dV Peak Tracking — 단계별 로드맵

> **2026-09 정확성 리뷰 추가:** 기존 peak 학습·추적 계획은 유지한다. [과학적 리뷰](SCIENTIFIC_REVIEW_2026-09.md)의 noise/smoothing artifact, phase 해석 제한, cell·batch hold-out 및 train-only golden/tuning 규칙을 적용한다. peak matching 성공이나 합성 ground truth 재현은 실제 열화 원인 증명이 아니다.

**갱신:** 2026-08-06  
**관련:** [ROADMAP.md](ROADMAP.md) Phase 1~2, [FEATURES.md](FEATURES.md) Tier 2, [GOLDEN_CYCLES.md](GOLDEN_CYCLES.md), [IMPROVEMENT_ROADMAP.md §12.5](IMPROVEMENT_ROADMAP.md#125-diffcapanalyzer에서-가져올-것-ica-피크-descriptor) (DiffCapAnalyzer peak descriptor)

> 목표: **피크를 믿을 수 있게 만들고 → 사이클마다 추적하고 → (데이터 충분 시) ML로 확장**  
> ICA descriptor 참고: peak V / H / area / W / sign / cycle-delta ([DiffCapAnalyzer](https://github.com/nicolet5/DiffCapAnalyzer) — 엔진 교체 없이 스키마 정렬)  
> **Full-cell 스택:** F1–F5는 하프셀 없이 완성 — [IMPROVEMENT_ROADMAP §0](IMPROVEMENT_ROADMAP.md#0-full-cell-우선-스택-하프셀-없이)

---

## Phase 0 — 지금 당장 (1~2일)

### Step 1. Peak review PNG 검수
- `example/docs/peak_review/` 아래 export PNG 확인
- 체크: 눈에 보이는 peak와 일치? 스파이크 오검출? P1/P2/P3 전압 순서?

### Step 2. Cycle 번호 규칙 고정
- LGES SJ900 set4: **Cycle-001 ≈ raw TotalCycle 2** (TC=1은 초기 Rest만)
- 이후 분석·라벨·golden 모두 **raw TotalCycle** 기준으로 통일 (classification 없을 때)

### Step 3. 기준 채널 확정
- **M01Ch025**를 1차 기준 채널로 사용 (2026-07-09 합의)
- Ch022/024는 Ch025 파이프라인 검증 후 확장

---

## Phase 1 — Golden 기준 만들기 (3~5일)

### Step 4. Golden cycle 지정
- “정상·저노이즈” curve를 golden으로 등록 → [GOLDEN_CYCLES.md](GOLDEN_CYCLES.md)
- Golden마다 charge/discharge **P1~P3 전압(V)·세기(H)** 기록
- 용도: peak assign window, 이상 기준, 추적 기대값

### Step 5. Peak 검출 파라미터 튜닝
- `cyclediag/features/dqdv_peaks.py` → `DqdvPeakConfig`
- Golden cycle에서 3+3 peak가 안정적으로 나올 때까지 조정
- prominence ↑ → 가짜 peak 감소 / ↓ → 약한 peak 누락

### Step 6. 동일 조건으로 다른 채널·사이클 export
- Ch024/025, 다사이클(10/50/100…) PNG + CSV 축적
- `export_dqdv_peak_review.py`, `score_cycle_dqdv_quality.py` 사용

---

## Phase 2 — 데이터 축적 (1~2주)

### Step 7. 라벨링 스키마
```text
cell_id, raw_total_cycle, leg, peak_num, V, H, phase_label, is_valid, note
```
- `phase_label`: 추후 `anode_1`, `cathode_H1`, `noise` 등
- 지금은 `is_valid` (yes/no)만이라도 기록

### Step 8. 다사이클 batch export
- BOL / 중간 / 후반 cycle에서 charge·discharge 각각 peak table 생성
- 목적: **V(t), H(t) 시계열** 확보 (열화 추적의 입력)

### Step 9. StepEnd 지표와 연결
- dQ/dV peak shift (ΔV, ΔH) + 용량 fade + Rest V를 한 feature table로 merge
- 어떤 peak 변화가 SoHQ와 같이 움직이는지 확인 → feature 우선순위

---

## Phase 3 — Peak Assign (규칙 + ML, 1주)

> **목표:** 매 사이클 검출된 bump에 **P1/P2/P3 이름(peak_id)** 을 붙인다.  
> 규칙(전압 밴드)이 기본, golden에서 학습한 **ML이 보조·대체**한다.

### Step 10. Golden 전압 window + centroid
- Golden cycle band-labeled peak에서 `(leg, peak_id)` 별 **V, H median** → `learned_criteria.json`
- 전압 window: median ± σ (기본 ±0.03~0.06 V)

### Step 11. Hungarian + ML peak assign (구현됨)
1. **학습:** golden cycle band-label → centroid(V,H window) + RF(`V`,`H`,`H_abs`)
2. **Hungarian:** `find_dqdv_peaks` 후보 ↔ peak_id **1:1 최적 매칭** (V·H cost)
3. **RF 가중:** cost × (1 − rf_proba) — 학습된 규칙이 매칭에 반영
4. **Hybrid:** 밴드 assign 후 **누락 peak_id**만 Hungarian으로 채움
5. **다셀 확장:** `train_peak_assign_multi()` + `tools/train_peak_assign_global.py`

| 모드 | 설명 |
|------|------|
| `band` | 전압 밴드만 |
| `hungarian` | Hungarian + learned centroid (ML 가중) |
| `ml` | `hungarian` 동일 |
| `hybrid` | **기본** — 밴드 + 누락분 Hungarian |

**구현:** `cyclediag/features/peak_assign.py` (`hungarian_assign_peaks`, `learned_criteria.json`)  
**글로벌 모델:** `example/docs/models/peak_assign_global_v1/`

### Step 12. Assign 품질 검증
- assign 실패율, peak missing rate, ΔV/ΔH 분포
- `assign_method=ml_*` 비율 모니터링
- 실패 구간만 수동 검수

---

## Phase 4 — 사이클 간 Peak 추적 (핵심)

> **세기(H)가 cycle마다 변해도 추적 가능**해야 함. 위치(V) 기준 assign + 세기는 별도 시계열로 관리.

### Step 13. 추적 단위 정의
- 입력 단위: `(cell_id, leg, peak_id)` × cycle
- `peak_id` = golden assign 결과 (P1/P2/P3 또는 phase_label)
- 추적량: `V_peak(cycle)`, `H_peak(cycle)`, optional `area_peak(cycle)`

### Step 14. 추적 방법 (권장 순서)

| 방법 | 원리 | 세기 변화에 강한 이유 |
|------|------|----------------------|
| **A. Voltage-anchor assign** | 매 cycle golden/이전 cycle V window로 peak 재매칭 | identity는 V에 묶음, H는 독립 추적 |
| **B. 정규화 세기** | `H_norm = H / H_golden` 또는 `H / ymax(cycle)` | 절대값 fade를 상대 비율로 |
| **C. Sequential prior** | cycle t 기대 V = V(t-1) + drift, window를 t-1 중심으로 이동 | golden만 고정하지 않고 궤적 추적 |
| **D. Local DTW / cross-correlation** | golden 주변 ±0.1 V slice에서 shape 매칭 | peak flatten 시에도 위치 보정 |
| **E. Kalman (선택)** | 상태 = (V, log H) per peak_id, cycle마다 update | 노이즈·missing peak 완화 |

### Step 15. 추적 feature table export
```text
cell_id, cycle, leg, peak_id, V, H, H_norm, dV_vs_golden, dH_vs_golden, assign_confidence
```
- **구현:** `cyclediag/features/peak_tracking.py` + `peak_export.export_peak_feature_table()`
- 산출: `{cell_id}_peak_tracking.csv`, `_peak_golden_ref.csv`, `_peak_tracking_summary.csv`, `plots/`

### Step 16. 이상·열화 연결 (구현됨)
- `peak_stepemd_join.py` — StepEnd SoHQ·용량 merge
- 산출: `{cell_id}_peak_cycle_merged.csv`, `_peak_fade_correlation.csv`
- 단일 셀: `dV/dcycle`, `dH/dcycle` 추세
- 다셀: 동일 cycle에서 peak feature 분산

---

## Phase 5 — ML 도입 (데이터 충분 후)

### Step 17. ML 도입 기준 (2개 이상 충족 시)
- [ ] 라벨 peak 3,000+
- [ ] 화학/공정/셀 타입 다양
- [ ] 규칙 assign 실패율 10%+
- [ ] batch마다 peak shape 크게 다름

### Step 18. 1차 ML = peak classifier (CNN 전보다 우선)
- Feature: V, H, width, leg, cycle, ΔV_golden, ΔH_golden
- Label: `phase_label`
- 모델: XGBoost / Random Forest
- **GroupKFold by cell_id**

### Step 19. 진단까지 연결
- 추적된 peak trajectory → `cyclediag` Phase 2 이상 탐지·분류

---

## 도구 매핑

| 단계 | 도구 / 모듈 |
|------|-------------|
| PNG 검수 export | `cyclediag/tools/export_dqdv_peak_review.py` |
| 저노이즈 cycle 추천 | `cyclediag/tools/score_cycle_dqdv_quality.py` |
| Golden 정의 | `planning/GOLDEN_CYCLES.md` |
| Peak 검출 | `cyclediag/features/dqdv_peaks.py` |
| Peak assign | `cyclediag/features/peak_assign.py` (Hungarian + RF, hybrid, multi-cell) |
| StepEnd join | `cyclediag/features/peak_stepemd_join.py` |
| Global assign train | `cyclediag/tools/train_peak_assign_global.py` |
| Trajectory export | `cyclediag/features/peak_export.py` (+ `peak_tracking.py`, `peak_plots.py`) |

---

## 한 줄 요약 (이번 주)

1. Ch025 golden 확정 ([GOLDEN_CYCLES.md](GOLDEN_CYCLES.md))
2. 전 cycle peak trajectory CSV 생성
3. V-anchor assign + `H_norm` 추적으로 세기 변화 추적
4. PNG는 이상 cycle만 spot-check
