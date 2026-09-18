# VP Diagnosis — Roadmap

**갱신:** 2026-08-06 · 버전: [VERSIONS.md](VERSIONS.md)

> **역사적 계획:** 아래 Phase 표는 초기 계획과 후속 연구 이력을 보존한다. 현재 코드는 추출·진단·보고서 및 회귀 테스트를 포함하므로 “최소 골격만 존재”하는 상태가 아니다. 구현과 과학 검증 상태는 [2026-09 리뷰](SCIENTIFIC_REVIEW_2026-09.md)를 우선 확인한다.

---

## 비전

PNE 사이클러 **Voltage Profile(VP)** 을 입력으로, 사이클·셀·공정 이상을 **자동 진단·스코어링**한다.

- pne_studio: **시각화·정량 분석** (dQ/dV, VP plot, mass, export)
- **cyclediag**: **패턴 인식·이상 탐지·분류** (ML + 통계 feature)

두 도구는 **데이터 형식만 호환**하고 코드베이스는 분리한다.

**성능·메모리 (합의 2026-07-01):** 다사이클 처리 시 pne_studio Calculate와 cyclediag 배치 extract 모두 **적응형 병렬** 정책 공유 — [§ 배치 성능](#배치-성능--병렬--메모리-pne_studio-연동).

---

## 추천 전략 (합의 전 기본안)

| 항목 | 추천 | 이유 |
|------|------|------|
| **1차 목표** | **이상 탐지 + golden reference 비교** | 라벨 부족해도 시작 가능 |
| **Feature** | **Tier 1·2 hand-crafted** 먼저 | 해석 가능, 데이터 적을 때 강함 |
| **모델** | Isolation Forest / One-Class SVM → XGBoost | 소규모 데이터·빠른 iteration |
| **딥러닝** | Phase 3 이후 (데이터 500+ 라벨) | VP 시퀀스 1D-CNN·DTW는 후순위 |
| **입력 단위** | **사이클 1개 × leg(charge/discharge)** | 라벨·해석 단위가 명확 |
| **UI** | CLI/배치 먼저 → Tk GUI는 Phase 4 | pne_studio와 역할 분리 |

---

## Phase 개요

```
Phase 0 ──► Phase 1 ──► Phase 2 ──► Phase 3 ──► Phase 4
 데이터·      Feature      Classical    Deep /       GUI·
 라벨 스키마   추출 MVP     ML           시계열 ML     리포트
```

### Degradation-mode 진단 로드맵 (별도 축)

상세: [LLI_LAM_DIAGNOSIS.md](LLI_LAM_DIAGNOSIS.md)

```
DM-P1 Full-cell pattern/est ──► DM-P2 온도 ──► DM-P3 Half-cell calibrate ──► DM-P4 통합
```

| DM Phase | 내용 | 상태 |
|----------|------|------|
| **1** | `vp_lges_cycle_v2` · LLI/LAM pattern score · confidence · calibration schema | **현재 범위** |
| **2** | 온도 동특성 · R 정규화 · 진단 안정성 | Planned |
| **3** | Half-cell OCV library · `*_est_hc_calibrated` · peak mapping · **DMA fit** (PyDMA 과학 + PyProBE API, §12.3–12.4) | Planned |
| **4** | full-cell + T + half-cell + post-mortem 통합 | Planned |

**원칙:** 하프셀 부재 ≠ LLI·LAM 진단 비활성. Level 1 pattern을 먼저 내고, Level 3은 검증·교정 인터페이스로 병행한다.

### Full-cell 우선 스택 (하프셀 없이 — 메인 트랙)

상세: [IMPROVEMENT_ROADMAP.md §0](IMPROVEMENT_ROADMAP.md#0-full-cell-우선-스택-하프셀-없이) · 실행 [§9.0](IMPROVEMENT_ROADMAP.md#90-full-cell-우선-하프셀-불필요--현재-메인-트랙)

```
F1 ICA/DVA → F2 detect → F3 match → F4 ΔV · F5 Δarea
         → F6 curve corr → F7 R/polarization → F8 change-point
         → pattern_score (fullcell_v1)
```

DM-P3 half-cell DMA는 위 스택 **완성 후** 교정용.

---

## Phase 0 — 데이터·라벨·평가 기준 (2주)

**목표:** 무엇을 “정상/이상”으로 볼지 정의하고, 샘플 세트를 모은다.

| ID | 작업 | 상태 |
|----|------|------|
| 0.1 | PNE CSV 컬럼 매핑 확정 | Planned |
| 0.2 | `DATA_SCHEMA.md` 샘플 파일 1~3개 수집 | Planned |
| 0.3 | Golden reference VP 정의 (정상 lot / baseline cycle) | Planned |
| 0.4 | 라벨 taxonomy 초안 (`LABELS.md`) | Draft |
| 0.5 | Train/val/test 분할 규칙 (셀 단위 hold-out) | Planned |

**산출물:** annotated CSV 또는 sidecar JSON 라벨 파일, 평가 체크리스트

---

## Phase 1 — Feature Engineering MVP (3~4주)

**목표:** VP 한 사이클 → 고정 길이 feature vector. **ML 없이**도 z-score·거리로 이상 후보 표시.

| ID | 작업 | 상태 |
|----|------|------|
| 1.1 | `io/cycler_csv.py` — PNE CSV 로드·충/방전 분리 | Scaffold |
| 1.2 | CC/CV 구간 분리 (전류 threshold) | Planned |
| 1.3 | Tier 1 features (용량, 평균 V, plateau, CE proxy) | Planned |
| 1.4 | Tier 2 features (dQ/dV peak, hysteresis, IR drop proxy) | Planned |
| 1.5 | Reference 대비 Mahalanobis / DTW distance | Planned |
| 1.6 | `features/extract.py` + 단위 테스트 | Scaffold |
| 1.7 | Feature parquet/CSV export CLI | Planned |

**산출물:** `python -m cyclediag extract --input ... --out features.parquet`

**Peak 추적 (병행):** [PEAK_TRACKING_ROADMAP.md](PEAK_TRACKING_ROADMAP.md) · [GOLDEN_CYCLES.md](GOLDEN_CYCLES.md)  
**열화·수명 확장:** [IMPROVEMENT_ROADMAP.md](IMPROVEMENT_ROADMAP.md) (§5.7 ΔQ(V), **§12** BatteryML·PyBaMM·**PyDMA·PyProBE**·DiffCapAnalyzer)

- `score_cycle_dqdv_quality.py` — 저노이즈 golden cycle 추천  
- (예정) `track_dqdv_peaks.py` — 전 cycle V/H trajectory

**배치 extract 시:** 사이클·파일 수가 많을 때만 병렬 (아래 [배치 성능](#배치-성능--병렬--메모리-pne_studio-연동)); 소수 사이클은 순차가 기본.

---

## Phase 2 — Classical ML (3~4주)

**목표:** 라벨이 있으면 분류, 없으면 이상 탐지. 해석 가능한 리포트.

| ID | 작업 | 상태 |
|----|------|------|
| 2.1 | Unsupervised: Isolation Forest, PCA + Hotelling T² | Planned |
| 2.2 | Supervised: XGBoost / Random Forest (라벨 있을 때) | Planned |
| 2.3 | SHAP 또는 permutation importance | Planned |
| 2.4 | Cross-validation (GroupKFold by cell_id) | Planned |
| 2.5 | `models/train.py`, `models/predict.py` CLI | Planned |
| 2.6 | 모델·scaler·feature spec 버전 고정 (joblib/pickle) | Planned |

**산출물:** `train` / `predict` / `evaluate` 서브커맨드

---

## Phase 3 — 시계열·딥러닝 (선택, 데이터 충분 시)

| ID | 작업 | 상태 |
|----|------|------|
| 3.1 | VP 리샘플링 (고정 N points, Q 또는 time 정규화) | Planned |
| 3.2 | 1D-CNN / LSTM autoencoder (이상 점수 = recon error) | Planned |
| 3.3 | Siamese network (reference pair distance) | Planned |
| 3.4 | 멀티사이클 입력 (cycle 1~N trajectory) | Planned |

---

## Phase 4 — 제품화 (4주+)

| ID | 작업 | 상태 |
|----|------|------|
| 4.1 | 배치 폴더 스캔 + HTML/PDF 리포트 | Planned |
| 4.2 | Tkinter GUI (roughness simulation / pne_studio 패턴) | Planned |
| 4.3 | pne_studio 연동 **없음** — 동일 CSV만 드래그앤드롭 | Planned |
| 4.4 | 모델 레지스트리 (버전·학습 데이터 해시) | Planned |

---

## 외부 오픈소스 참고 (BatteryML · PyBaMM · PyDMA · PyProBE · DiffCapAnalyzer)

상세: [IMPROVEMENT_ROADMAP.md §12](IMPROVEMENT_ROADMAP.md#12-외부-오픈소스-참고-batteryml--pybamm)

**원칙:** 패키지 통째 내장 금지. 과학 기능·API 계약·feature·평가만 흡수.  
**DMA 추천:** PyDMA(과학) + PyProBE(API 구조).

| 출처 | 가져올 아이디어 | cyclediag 반영 위치 |
|------|-----------------|---------------------|
| [BatteryML](https://github.com/microsoft/BatteryML) | Cycler별 파서 → `CellData`/`CycleData`/`FeatureSet` | `io/adapters/`, schema |
| BatteryML | YAML 실험 · feature↔model · train/val/test | Phase 2 evaluate |
| BatteryML | Severson ΔQ(V) | §5.7, RUL |
| [PyBaMM](https://github.com/pybamm-team/PyBaMM) | 전극 SOH 출력 계약 | Layer 2 schema |
| [PyDMA](https://github.com/tum-ees/PyDMA) | LLI/LAM_NE/LAM_PE · stoich window · blend · OCV+DVA+ICA 가중 fit | DM-P3, `diagnosis/halfcell/` |
| [PyProBE](https://github.com/ImperialCollegeLondon/PyProBE) | `OCP`/`CompositeOCP` · fit target · batch DMA · warm-start · `quantify_degradation_modes` | 동일 모듈 API |
| [DiffCapAnalyzer](https://github.com/nicolet5/DiffCapAnalyzer) | ICA peak V/H/area/W/sign · cycle delta · baseline | peak schema, review UI |

**Phase 연계 (초안):**
- Phase 1+: FeatureSet + ΔQ(V) + ICA peak descriptor 정렬
- Phase 2+: YAML experiment + RUL
- DM-P3: OCP/CompositeOCP + batch DMA (하프셀·템플릿 확보 후 수치 채움)
- PyBaMM: 선택적 forward 합성만 (런타임 비의존)

---

## Feature Tier 요약

상세: [FEATURES.md](FEATURES.md)

| Tier | 예시 | Phase |
|------|------|-------|
| **0** | Raw Q, V, t, I 시퀀스 | 3 |
| **1** | Q_max, V_mean, CC fraction, CV time | 1 |
| **2** | dQ/dV peak V/H/W, ΔV hysteresis, dV/dt | 1–2 |
| **3** | Embedding (CNN), DTW to golden | 3 |

---

## 리스크·완화

| 리스크 | 완화 |
|--------|------|
| 라벨 부족 | Phase 1 reference distance + unsupervised 먼저 |
| 전극·공정별 VP 형태 차이 | `metadata.cathode`, `loading` 으로 그룹별 reference |
| CV 구간 노이즈 | CC-only feature 옵션 + smoothing 파라미터 고정 |
| pne_studio 중복 | cyclediag는 진단만; plot은 export PNG를 리포트에 첨부 |
| 소수 사이클에서 병렬 | 워커 기동·직렬화 오버헤드 → **적응형** (임계값 미만은 순차) |
| 대량 사이클 RAM | pne_studio Memory Phase(v1.10) 선행 · cyclediag는 feature만 유지 |

---

## 배치 성능 · 병렬 · 메모리 (pne_studio 연동)

cyclediag는 **사이클×leg 단위 feature 추출**을 배치로 돌린다. pne_studio **Calculate & Plot**과 같은 CSV·사이클 단위이므로, **병렬·메모리 정책을 맞춘다.**

### 병렬 계산 — 적응형 (항상 ON 아님)

| 사이클 수 (대략) | 병렬 이득 | 권장 |
|------------------|-----------|------|
| **1~15** | 거의 없음 · **더 느려질 수 있음** | **순차** |
| **15~40** | 파일·포인트 밀도에 따라 애매 | 순차 또는 자동 |
| **40~100+** | CPU 멀티코어 활용 가능 | 병렬 (워커 제한) |
| **200+** | RAM 피크가 병목 | 메모리 최적화 + 워커 2~4 |

**이유:** 프로세스/풀 기동, 입력 pickle, 결과 수집 등 **고정 오버헤드**가 사이클당 계산 시간보다 크면 소수 N에서 손해.

**구현 원칙 (pne_studio Calculate · cyclediag `extract` 공통):**

| 항목 | 값 |
|------|-----|
| 기본 | **순차** (현재 pne_studio 동작과 동일) |
| 자동 전환 | `N_cycles >= 임계값` (초안 **20~30**, 벤치 후 확정) |
| 워커 수 | `min(cpu_count - 1, 4)` — RAM 피크가 코어 수만큼 늘지 않게 |
| UI/CLI | 옵션: `Parallel` · `Auto` (default) · `Off` |
| Phase | pne_studio **v1.10+** (Memory M1–M3 후) · cyclediag **Phase 1.7+** |

> **cyclediag Phase 1** 본체는 ML 없이 feature만 — 병렬은 **추출 속도**용. **ML 학습·추론(Phase 2+)** 은 별도 배치이나 동일 적응형 원칙 적용 가능.

### 메모리 — pne_studio Memory Phase (v1.10)

병렬만으로는 **RAM 부족**을 해결하지 못함. 대량 사이클 전에 pne_studio 쪽 최적화를 권장.

| ID | 내용 | cyclediag와의 관계 |
|----|------|------------------|
| M1 | 로드 시 `float32` downcast | 동일 CSV 읽기 부담 ↓ |
| M2 | VP 캐시 다운샘플 (`vp_orig` 제거) | Studio 캐시 축소 · cyclediag는 Tier1/2만 export |
| M3 | 계산 경로 `.copy()` 감소 | Calculate·extract 공통 이득 |
| M4 | Calculate 후 raw CSV unload (옵션) | 플롯은 캐시 · 재계산 시 reload |
| M5 | 25사이클마다 GC · 200+ 경고 | OOM 예방 |

**문서:** `pne_studio/planning/specs/memory-v2.md` · 릴리스 번호 **v1.10** ([pne_studio ROADMAP](../../pne_studio/planning/ROADMAP.md))

**cyclediag 배치:** feature parquet만 남기고 **원본 VP 시퀀스는 기본 저장 안 함** (Tier 0는 Phase 3). → Studio보다 메모리 footprint 작게 유지.

### 우선순위 (성능)

1. **pne_studio v1.10** Memory M1–M5 (RAM · copy)
2. **적응형 병렬** Calculate / `cyclediag extract` (사이클 많을 때만)
3. Phase 2+ lazy load · Parquet cache (pne_studio Memory Phase 2, cyclediag 대용량 배치)

---

## 의사결정 backlog (NOTES.md로 이동)

- [x] 병렬 Calculate/extract — **적응형**, 소수 사이클은 순차 (2026-07-01)
- [ ] 병렬 임계값 확정 (20 vs 30) — 실측 벤치 후
- [ ] 첫 타깃: **공정 QC** vs **수명 예측** vs **불량 분류**
- [ ] Charge only vs charge+discharge feature
- [ ] 사이클 번호: formation만 vs 전체 life
- [ ] 외부 라이브러리: `scikit-learn`, `xgboost`, `pyod` (이상탐지)

## 2026-09-18 Review — 과학적 정확성과 실행 순서

[과학적 계약 및 외부 검증 계획](SCIENTIFIC_REVIEW_2026-09.md)에 이번 리뷰의 구현됨/계획/실험적 상태와 단위·CE/VE·IR·peak·품질 게이트 한계를 정리했다. 기존 상세 알고리즘 제안은 삭제하지 않는다.

우선순위는 **데이터 계약 → 검증된 feature 추출 → 독립 측정으로 외부 검증한 진단 → 재현 가능한 batch/cloud·보고서**다. Full-cell proxy와 미교정 confidence를 원인 확률로 승격하지 않는다. 데이터 부족은 unknown으로 보존하며, cell 및 batch hold-out·독립 용량/OCV/조건에 맞는 half-cell/EIS 검증과 result hash/config 영속화는 별도 인수 대상이다.
