# cyclediag 개선 로드맵 · 통합판

> **2026-09-18 reproducibility update:** bounded input/config/code provenance and checked feature CSV sidecars are implemented for selected high-level API/CLI entry points. See [implemented scope and remaining work](REPRODUCIBILITY.md); older all-pipeline manifest proposals below are not claims of complete coverage.

> **2026-09 정확성 리뷰 추가:** 아래 연구 제안과 이력은 보존한다. 현재 단위·품질 게이트의 실제 범위, full-cell/ECM 비식별성, 단계별 독립 검증 기준은 [과학적 리뷰](SCIENTIFIC_REVIEW_2026-09.md)를 함께 적용한다. “구현됨”은 메커니즘 식별이나 외부 검증 완료가 아니다.

> **대상 셀:** ASSB SJ900 (S83S 양극 / SJ-ASG903-1300 Si-rich 음극), 2.5–4.2 V, 45 °C
> **프로토콜:** routine 0.5C CC-CV · C/3 RPT 2사이클 (약 105 사이클 주기) · DC-IR (SOC 20/50/80, 방전, 1C, 30 s)
> **기준 데이터:** SJ900 set4 Ch22 (564 cycles, SoHQ ~65 %), Ch25
> **갱신:** 2026-08-06 — **§0 Full-cell 우선 스택** 추가 (하프셀 없이 ICA/DVA·peak·R·change-point). §12 OSS 참고 유지

---

## 목차

- [0. Full-cell 우선 스택 (하프셀 없이)](#0-full-cell-우선-스택-하프셀-없이)
- [1. 확정 사실 요약](#1-확정-사실-요약)
- [2. 구조적 한계 진단](#2-구조적-한계-진단)
- [3. ASSB 열화 모드 체계 재정의](#3-assb-열화-모드-체계-재정의)
- [4. 지표 체계 정비](#4-지표-체계-정비)
- [5. 알고리즘 명세](#5-알고리즘-명세)
  - [5.1 dQ/dV 필터 파라미터 스윕 진단](#51-dqdv-필터-파라미터-스윕-진단)
  - [5.2 적응형 스무딩 (V-폭 기준)](#52-적응형-스무딩-v-폭-기준)
  - [5.3 R(t) 3성분 분해](#53-rt-3성분-분해)
  - [5.4 자가방전율 / 미세단락 검출](#54-자가방전율--미세단락-검출)
  - [5.5 펄스 후 회복 τ](#55-펄스-후-회복-τ)
  - [5.6 3-파라미터 곡선 정합](#56-3-파라미터-곡선-정합)
  - [5.7 ΔQ(V) 통계](#57-δqv-통계)
  - [5.8 Q-domain 제약 디컨볼루션 + 피크 그룹](#58-q-domain-제약-디컨볼루션--피크-그룹)
  - [5.9 η(SOC) × DC-IR 결합 R_eff](#59-ηsoc--dc-ir-결합-r_eff)
  - [5.10 Q_relax · RCF · PER](#510-q_relax--rcf--per)
  - [5.11 SOC 분해 히스테리시스](#511-soc-분해-히스테리시스)
  - [5.12 페이드 지수 · knee 검출](#512-페이드-지수--knee-검출)
  - [5.13 데이터 품질 게이팅](#513-데이터-품질-게이팅)
  - [5.14 CC/CV 신호 기반 재검출](#514-ccv-신호-기반-재검출)
- [6. 파일 IO 견고화](#6-파일-io-견고화)
- [7. 아키텍처 개선](#7-아키텍처-개선)
- [8. 검증 체계](#8-검증-체계)
- [9. 실행 계획](#9-실행-계획)
- [10. 미해결 질문](#10-미해결-질문)
- [11. 금지 사항](#11-금지-사항)
- [12. 외부 오픈소스 참고 (BatteryML · PyBaMM · PyDMA · PyProBE · DiffCapAnalyzer)](#12-외부-오픈소스-참고-batteryml--pybamm)


---

## 0. Full-cell 우선 스택 (하프셀 없이)

> **제품 결정 (2026-08-06):** 하프셀 OCP가 없어도 **진단·추적을 최대한 완성**한다.
> Half-cell DMA(PyDMA/PyProBE, §12.3–12.4)는 **검증·교정 레이어**로만 뒤에 둔다.
> 출력은 `diagnosis_version=fullcell_v1` + `*_pattern_score` / proxy estimate.
> **`*_est_hc_calibrated`는 하프셀·검증 템플릿 전까지 채우지 않는다.**

### 0.1 지금 해야 할 8가지 (핵심)

| # | 기능 | 하프셀? | 상태 (대략) | 주 모듈 / § | 진단에 쓰는 방식 |
|---|---|---|---|---|---|
| F1 | **ICA/DVA 곡선 생성** | 불필요 | 부분 구현 | `dqdv_peaks` · `dqdv_segment` · §5.1–5.2 | ICA=`dQ/dV`, DVA=`dV/dQ` 둘 다 export·캐시 |
| F2 | **피크 자동 검출** | 불필요 | 구현됨 | `dqdv_peaks.py` · DiffCapAnalyzer descriptor §12.5 | V, H, area, W, sign |
| F3 | **피크 matching** | 불필요 | 부분 구현 | `peak_assign` · `peak_tracking` · `peak_evolution` | golden/RPT 앵커 · Hungarian / Viterbi |
| F4 | **피크 위치 이동** | 불필요 | 부분 구현 | `peak_trajectory` · evolution | ΔV, ΔQ → LLI / slippage **pattern** |
| F5 | **피크 면적 감소** | 불필요 | 부분 구현 | group_area · peak area trajectory · §5.8 | → LAM_PE **pattern** (Si: LAM_NE는 피크 단독 금지) |
| F6 | **곡선 correlation** | 불필요 | 일부 | ΔQ(V) §5.7 · DTW/corr · fade_correlation | baseline vs cycle-N 형상 유사도 |
| F7 | **저항·polarization 증가** | 불필요 | 부분~계획 | DC-IR · §5.3–5.5 · §5.9–5.10 PER | 접촉/계면/확산 분리 (ASSB 핵심) |
| F8 | **change-point 탐지** | 불필요 | 계획 | §5.12 knee · §7.7 PELT/Bacon-Watts | fade·R·peak 궤적의 급변 시점 |

```
full-cell raw
  → F1 ICA/DVA
  → F2 detect → F3 match → F4 ΔV · F5 Δarea
  → F6 curve corr / ΔQ(V)
  → F7 R · polarization
  → F8 change-point / knee
  → pattern_scoring (LLI / LAM_PE / impedance / contact …)
  → (later) half-cell DMA calibrate   ← optional, not blocking
```

### 0.2 하프셀 없이 가능한 것 / 보류할 것

| 가능 (지금) | 보류 (하프셀·OCP 템플릿 후) |
|---|---|
| ICA/DVA · peak detect/match/track | Aged half-cell 대비 `*_est_hc_calibrated` 검증 |
| peak ΔV / Δarea → **pattern score** | stoichiometry window 절대값 (열화 후 OCP 필요) |
| ΔQ(V) · curve corr → early fade / RUL proxy | blend phase OCP 분해 (aged) |
| R 3성분 · PER · hysteresis | electrode utilization 절대 % (aged) |
| knee / change-point | full-cell↔half-cell peak 화학 라벨 확정 |
| 3-param curve fit proxy (§5.6) | |
| **BOL 하프셀 OCP library** (`example/fixtures/halfcell/`, C/20) | **Aged / harvested 하프셀 — 아직 없음** |

### 0.3 완료 정의 (Full-cell MVP DoD)

- [ ] 사이클마다 ICA **및** DVA 곡선 parquet/npz 저장 (동일 전처리 버전 해시)
- [ ] 피크 테이블: match_id, V, H, area, W, sign, confidence
- [ ] 궤적: `dV_vs_baseline`, `dArea_vs_baseline`, match_confidence
- [ ] 곡선: `corr_to_baseline`, `dQV_log_var` (또는 equiv)
- [ ] 저항: `R_ohmic`/`R_ct`/`A_diff` 또는 가용 DC-IR proxy + PER
- [ ] change-point: SoHQ·R·대표 peak_V 중 ≥1 시계열에 변점 + knee flag
- [ ] `*_pattern_score`가 위 증거를 `supporting_features`로 인용
- [ ] 하프셀 없이도 CLI 한 방으로 셀 리포트 생성

### 0.4 실행 순서 (이 스택만 뽑은 것)

```
1. F1+F2 품질 고정   (§5.1 스윕 → §5.2 스무딩)     ← 피크 전부이 선행
2. F3 matching 안정  (RPT 앵커 · evolution)
3. F4+F5 궤적 export (ΔV, Δarea → diagnosis 입력)
4. F7 R/polarization (§5.3 우선 — ASSB 접촉 손실)
5. F6 curve corr     (§5.7 ΔQ(V) + baseline corr)
6. F8 change-point   (§5.12 + §7.7)
7. pattern_scoring 재배선 (mode_weights_assb_si_v1)
```

상세 알고리즘은 §5·§7, peak 단계는 [PEAK_TRACKING_ROADMAP.md](PEAK_TRACKING_ROADMAP.md),
정책은 [LLI_LAM_DIAGNOSIS.md](LLI_LAM_DIAGNOSIS.md).


## 1. 확정 사실 요약

`IMPROVEMENT_ANSWERS.md`로 확인된 내용 중 설계에 직접 영향을 주는 항목.

### 1.1 유리한 조건 (즉시 활용 가능)

| 항목 | 값 | 열리는 기능 |
|---|---|---|
| DC-IR 샘플링 | **전 구간 0.1 s (10 Hz)**, t≤1 s 11점 | R 3성분 분해 (§5.3) |
| DC-IR raw 저장 | **raw V–I 트레이스 보존** | 동일 |
| 펄스 전 rest | **3600 s × 3 SOC** | 자가방전율 (§5.4) |
| 펄스 후 rest | **1800 s** | R_recovery_τ (§5.5) |
| RPT 2사이클 | Δ = 0.045/68.9 Ah ≈ **0.065 %** | Q_relax 노이즈 하한 확정 (§5.10) |
| `POST_RPT_EXCLUDE` | **5, 이미 구현** | 회복 구간 오염 방지 완료 |
| Ch22 수명 | **564 cycles, SoHQ ~65 %** | knee 검증 데이터 (§5.12) |
| 두 rate 곡선 | C/3 + 0.5C 동일 셀 | RCF, η(SOC), PER (§5.9–5.10) |
| **하프셀 BOL OCP** | **있음** — `example/fixtures/halfcell/` (C/20, 음극 cycle 1–3) | DM-P3 OCP library · DMA fit 프로토타입 |

### 1.2 제약 조건 (설계에 반영 필요)

| 항목 | 값 | 영향 |
|---|---|---|
| **화학** | ASSB, Si-rich 음극 | 모드 체계 전면 재정의 (§3) |
| **온도 로그** | `Temp (Celsius)` 전부 0.0 | Arrhenius 보정·DTV 불가. 결론은 45 °C 한정 |
| **하프셀 aged OCP** | **없음** (BOL만 있음) | 열화 전후 OCP 비교·harvest 검증 불가. BOL로 library·fit 프로토타입은 가능 |
| **golden set** | 없음 | 합성 데이터 + 공개셋 검증이 유일 (§8) |
| **RPT 주기** | 105 사이클 | 앵커 간격 과다 → §5.8 보간 신뢰도 리스크 |
| **DC-IR SOC 배치** | 3개 **연속 사이클**에 분산 | SOC 간 비교 시 순서 효과 존재 |
| **펄스 전류** | 1C (77.34 A) @ 45 °C | 선형 응답 이탈 가능 → 절대값 대신 상대 추이 |
| **SOC 기준** | 직전 capa 대비 **상대 %** | 열화 시 측정점 이동이 R 추이에 혼입 |
| **구속 압력** | UNKNOWN | ASSB 접촉 손실 해석의 지배 변수 (최우선 확인) |
| `chgCVcapa` | raw에 있으나 **feature=0** | CV 지표 전체 무효 → 버그 (§5.14) |

### 1.3 현행 dQ/dV 설정을 물리 단위로 환산

| 설정 | 값 | 환산 (72 Ah 기준) |
|---|---|---|
| `n_interp` | 500 | **0.144 Ah / 점** |
| `interp_axis` | `"Q"` | 평탄 구간에서 V 해상도 과밀 |
| `sg_window` | 21 (도구 31) | **≈ 3.0 Ah** (31이면 4.5 Ah) |
| `min_distance_frac` | 0.04 | 20점 = **≈ 2.9 Ah** |
| merge 임계 | 12 mV | — |

평탄 구간에서 3 Ah 폭 스무딩은 전압으로 **수십 mV** 규모다.
**스무딩 폭이 병합 임계값과 같은 자릿수이거나 크다** → 0.5C 피크 병합이
물리 현상이 아니라 처리 인공물일 가능성이 높다. §5.1이 이를 판정한다.

---

## 2. 구조적 한계 진단

### 2.1 역문제를 풀지 않는다

`pattern_scoring.py`는 **증거 → tanh → 가중합 → 점수**의 단방향 매핑이다.

```
현재:  feature 값들  --(가중치)-->  LLI_pattern_score = 0.72
필요:  측정 V(Q)     --(fitting)--> 파라미터 추정 + 공분산
```

역문제로 전환하면 heuristic "mode collision 감점"이 **파라미터 공분산 행렬**로 대체된다.
두 모드가 구분 안 되면 신뢰도를 깎는 게 아니라 **원리적 분리 불가**를 수학적으로 출력한다.

§5.6(곡선 정합)이 반쪽셀 없이 가능한 1차 근사이고, 반쪽셀 확보 후 완전한 역문제로 확장한다.

### 2.2 검증 수단이 없다

golden set이 없어 heuristic 튜닝의 과적합을 판별할 수 없다.
→ §8의 합성 데이터 역복원이 유일한 즉시 가능한 검증 경로다.

### 2.3 지표 폭발 (188 columns)

동일 물리량이 raw / delta / slope / ratio / mode evidence / anomaly input으로 중복 계상된다.
z-score anomaly가 상관된 지표를 중복 카운트해 신호를 왜곡한다.
→ §4.1의 패밀리 태그로 anomaly 입력을 제한한다.

### 2.4 품질 점수가 실측 근거가 없다

현행:
```
mode data_quality = n_available_evidence / n_configured_terms
row score = mean(mode data_quality)
confidence = 0.45*dq + 0.40*agree + 0.15*min(1, n/5)
```
이는 **"증거 컬럼이 채워졌는가"** 만 본다. 노이즈, 샘플 밀도, 피팅 잔차 같은
실제 데이터 품질은 전혀 반영하지 않는다. → §5.13으로 교체.

---

## 3. ASSB 열화 모드 체계 재정의

### 3.1 무효화 항목 (`mode_weights`에서 제거)

| 항목 | 무효 사유 |
|---|---|
| 전해액 고갈 / wetting 저하 | **액체 전해질 없음** |
| 액상 확산 제한 | 동일 |
| `Q_NE_apparent` (흑연 stage 간격) | **Si는 stage 전이 피크 없음.** 비정질화 후 특징 소실 |
| 피크 기반 LAM_NE 분해 | 관측 피크는 **대부분 양극 유래** |
| 충전 vs 방전 곡선 직접 비교 | Si 히스테리시스, **OCV 경로 의존** |
| 통상적 Li plating 신호 | ASSB에서는 덴드라이트 관통으로 발현, 신호 형태가 다름 |

### 3.2 재정렬된 모드 우선순위

| 순위 | 모드 | ASSB에서의 의미 | 주 증거 | 산출 §
|---|---|---|---|---|
| 1 | **계면 저항 증가** (SE\|CAM, SE\|AAM) | ASSB 1차 열화 경로 | `R_ct`, `tau_ct` 성장 | 5.3 |
| 2 | **화학-기계적 접촉 손실** | Si 팽창·수축 → 계면 박리, void | **`R_ohmic` 성장** | 5.3 |
| 3 | **SE 분해** | 계면 전위창 초과 (45 °C 가속) | CE 저하 + `R_ct` 동시 성장 | 5.3 + CE |
| 4 | **덴드라이트 관통 / 미세단락** | SE 결정립계 침투 | **자가방전율** | 5.4 |
| 5 | LAM_PE | 양극 활물질 손실 | 피크 면적·간격 (양극 유래이므로 유효) | 5.8 |
| 6 | LLI | 리튬 재고 손실 | CE, 슬리피지, 곡선 offset | 5.6 |
| 7 | 고상 확산 제한 | SE 내부 + 입자 내 | `A_diff`, PER | 5.3, 5.10 |

> **`contact_loss_score`를 최상위로 이동한다.**
> ASSB에서 이것은 예외적 실패가 아니라 **주 열화 경로**다.

### 3.3 R 성분의 물리적 귀속 — 액체 셀과 다름 (핵심)

| 성분 | 액체 Li-ion | **ASSB (본 셀)** |
|---|---|---|
| `R_ohmic` ↑ | 전해액 열화 (느림, 부차적) | **유효 접촉 면적 손실 · void → 접촉 손실 직접 지표** |
| `R_ct` ↑ | SEI 성장 | 계면상 성장 · SE 분해 |
| `A_diff` ↑ | 액상 + 고상 확산 | **고상 확산만** (SE + 입자 내) |
| `R_ohmic / R_total` | 대체로 안정 | **상승 = 기계적 열화 지배 신호** |

ASSB에서 R 3성분 분해는 부가 기능이 아니라
**접촉 손실과 계면 화학 열화를 구분하는 유일한 저비용 수단**이다.

### 3.4 3계층 온톨로지 (출력 스키마에 반영)

```
Layer 1  관측 (observation)  : 용량 손실, 분극 증가, CE 저하, 피크 이동, 자가방전
   ↓  [곡선 정합 · R 분해로 비교적 확실하게 이동 가능]
Layer 2  상태 (state)        : LLI, LAM_PE, R_ohmic/R_ct/R_diff 증가
   ↓  [거의 항상 추가 실험 필요 — UI에 명시]
Layer 3  기구 (mechanism)    : 접촉 손실, SE 분해, 관통, 입자 균열
```

출력 컬럼에 `layer` 속성을 부여하고, UI는 Layer 3을 **가설**로 표기한다.

### 3.5 Si-rich 음극 대응

| 항목 | 조치 |
|---|---|
| SOC 분해 히스테리시스 | **최우선 지표로 격상** — Si 열화 직접 지표 (§5.11) |
| 곡선 정합 | 충·방전 **각각 자기 방향 baseline에만** 정합. 교차 비교 금지 |
| 피크 해석 | **LAM_PE 전용**. LAM_NE는 피크로 판정하지 않음 |
| NP ≈ 1.08 | Si-rich 치고 타이트 → 저 SOC 저항 급증 시 국부 과리튬화 맥락 |
| 팽창/수축 | 접촉 손실의 주 구동력 → `R_ohmic` 추이와 연결 해석 |

---

## 4. 지표 체계 정비

### 4.1 패밀리 태그 (선행 필수) — **구현됨** (§4.4)

지표를 늘리기 전에 정리한다. 상관된 지표를 더 넣으면 z-score anomaly가 악화된다.

```python
FAMILY = {
    "coulombic":     ["dchgCapa", "SoHQ", "CE", "CI", "k_SEI", ...],
    "kinetic":       ["R_ct", "tau_ct", "R_DCIR_*", "tau_CV", ...],
    "geometric":     ["R_ohmic", "R_ohmic_frac", "contact_*", ...],
    "transport":     ["A_diff", "PER", "eta_SOC*", ...],
    "thermodynamic": ["peak_*_V", "peak_*_Q", "fit_offset", "fit_scale", ...],
    "integrity":     ["self_discharge_rate", ...],
    "thermal":       ["dT_max", "DTV_*"],   # 온도 로그 복구 후
}
```

**규칙**
1. anomaly 입력은 **패밀리당 대표 1~2개**만
2. 나머지는 screening · diagnosis 전용
3. 상관행렬 산출해 |r| > 0.95 쌍은 빌드 시 경고
4. registry에 `unit`, `layer`, `family`, `requires`, `valid_when` 기록

위 초안은 물리 계열(coulombic/kinetic/…)로 묶는 안이었으나, 실측 결과 **계열이
아니라 개별 측정량 단위로 묶어야** 함이 확인됐다. 예를 들어 `hyst_area`(global)와
`hyst_area_low`는 같은 "coulombic/thermodynamic" 계열이지만 |r| = 0.004로 완전히
독립이다. 반대로 `hyst_area`와 `hyst_max_dV`는 같은 dV(Q) 구간의 두 요약값이라
per-cell 0.988로 중복이다. 따라서 실제 구현(§4.4)은 계열 대신 **측정량 패밀리**를
쓰고, 계열 정보는 `role`로 분리했다.

### 4.2 신규 지표 총괄

| 지표 | 정의 위치 | 모드 대응 | 필요 데이터 | 난이도 |
|---|---|---|---|---|
| `R_ohmic`, `R_ct`, `tau_ct`, `A_diff` | §5.3 | 접촉/계면/확산 분리 | 10 Hz 펄스 ✅ | 낮음 |
| `mech_vs_chem_ratio` | §5.3 | 기계적 vs 화학적 우세 | 위와 동일 | 낮음 |
| `self_discharge_rate(SOC)` | §5.4 | 관통·미세단락 | rest 3600 s ✅ | 매우 낮음 |
| `R_recovery_tau` | §5.5 | 확산 교차검증 | rest 1800 s ✅ | 낮음 |
| `fit_scale/offset/dR/residual` | §5.6 | LAM/LLI/R 분해 | routine 곡선 ✅ | 중간 |
| `dQV_var` 외 ΔQ(V) 통계 | §5.7 | 조기 수명 예측 | routine 곡선 ✅ | 낮음 |
| `group_area/centroid/width` | §5.8 | LAM_PE | dQ/dV ✅ | 높음 |
| `eta_SOC20/50/80`, `eta_argmax` | §5.9 | 제한 전극 식별 | RPT+routine ✅ | 낮음 |
| `R_ratio_20_50`, `R_SOC_slope` | §5.9 | SOC 의존 형상 | DC-IR 3점 ✅ | 매우 낮음 |
| `Q_relax` | §5.10 | 완화 성분 축적 | RPT 2사이클 ✅ | 매우 낮음 |
| `rate_capability_fade` | §5.10 | kinetic vs thermo | RPT+routine ✅ | 낮음 |
| `polarization_excess_ratio` | §5.10 | 확산 제한 | +DC-IR ✅ | 낮음 |
| `hyst_area_SOC*` | §5.11 | **Si 열화** | routine 곡선 ✅ | 낮음 |
| `fade_exponent_b`, `knee_*` | §5.12 | 지배 기구 판별 | SoHQ 궤적 ✅ | 중간 |
| `VE`, `EE` | §4.3 | 분극 vs 쿨롱 분리 | 에너지 적분 ✅ | 매우 낮음 |
| `CI_per_hour` | §4.3 | 캘린더 vs 사이클 | 시간 ✅ | 매우 낮음 |
| `samples_per_mV`, `v_noise_sigma` 등 | §5.13 | 품질 게이팅 | raw ✅ | 낮음 |
| `tau_CV`, `Q_CV_norm` | §5.14 | 음극 수용 한계 | CV 검출 수정 후 | 중간 |

### 4.3 즉시 추가 가능한 단순 지표

```python
# 에너지 효율 계열 — 분극 손실과 쿨롱 손실 분리
E_chg  = ∫ V·I dt  (충전 leg)
E_dchg = ∫ V·I dt  (방전 leg)
EE = E_dchg / E_chg
VE = EE / CE          # 전압 효율. 분극 손실만 반영
dE = E_chg - E_dchg   # 비가역 발열 프록시

# 쿨롱 비효율 시간 정규화 — 45 °C 캘린더 성분 분리
# CI = 100 - CE 자체는 컬럼으로 내보내지 않는다 (§4.4: CE의 아핀변환, 정보량 0).
# 시간 정규화된 형태만 내보낸다 — CE와 pooled 0.90 / per-cell 0.96 으로 독립.
CI_per_hour = (100 - CE) / cycle_duration_hours
CI_cumulative = cumsum(100 - CE)
k_SEI = fit(Q_loss ~ a·sqrt(N))   # a 계수
```

`VE`는 ASSB에서 특히 유용하다. 계면 저항 증가는 CE를 거의 건드리지 않고 VE만 떨어뜨리므로,
`CE 유지 + VE 하락` 조합이 **계면 열화의 깨끗한 신호**가 된다.

### 4.4 중복 지표 통폐합 — **구현됨**

모듈: `cyclediag/features/indicator_registry.py`

#### 문제

extract 파이프라인은 사이클당 393개 컬럼을 냈다. 그중 상당수가 같은 측정량의
다른 이름·단위·정규화였다. "숫자형 컬럼 전부"로 feature pool을 만드는 소비자
(anomaly score, indicator screen, regime singularity, lifetime ranking)는
**어떤 물리량이 별칭을 몇 개 가졌는지에 따라 그 물리량의 가중치가 결정**됐다.
실측된 최악 사례:

| 패밀리 | 구 pool에서 중복 계상 | 컬럼 |
|---|---|---|
| `rest_v_eoc_relaxed` | ×6 | `EoC_restV_{60s,30m,end}` + 각 `delta_` |
| `capacity_discharge` | ×4 | `dchgCapa`, `SoHQ`, `dchg_E`, `dchg_dVdQ_SOC0_Q` |
| `capacity_charge` | ×4 | `chgCapa`, `chgCCcapa`, `chg_E`, `chg_dQdV_area_sum` |
| `cv_amount` | ×4 | `chgCVcapa`, `chgCVtime`, `chgCapa_CCratio`, `delta_chgCapa_CCratio` |
| `r_eod_chg_10s` | ×4 | `EoD_chgR_10s`, `_T25`, `_inc`, `_R10_minus_R0p1` |

구 prefix 필터 pool 168개 컬럼 = 119개 패밀리, **32개 패밀리가 중복 계상**,
중복 슬롯 49개. 통폐합 후 167개 컬럼 = 167개 패밀리, **중복 0**.

#### 두 개의 직교축

- `role` — 양의 종류. `indicator`(열화 증거) / `target`(건강도 및 그로부터
  계산된 양) / `covariate`(프로토콜·환경) / `qc`(provenance) / `diagnosis`(진단
  출력) / `meta`. **자동 pool에는 `indicator`만 들어간다.**
- `family` — 물리 측정량. 패밀리당 대표 1개만 자동 pool에 진입한다. 나머지
  멤버는 계속 emit되고 명시적 config·플롯·export에서 그대로 쓸 수 있다.
  패밀리는 "버려라"가 아니라 **"독립으로 세지 말라"**는 뜻이다.

`delta_<x>` / `<x>_inc`는 셀 내에서 `<x>`의 아핀변환이므로 `canonical_base`가
자동으로 접어 넣는다 (패밀리 목록에 별도 등재하지 않음).

#### 통폐합 판정 기준

(1) 또는 (2)를 만족하고, **항상 (3)을 만족**해야 병합한다.

1. **구조적 동일성** — 한쪽이 다른 쪽의 정확한 대수 함수. 단위 변환
   (`f_Q_spec` = `f_Q_max` / active mass), 같은 양을 느리게 변하는 기준으로
   정규화한 형태 (`Q_relax_pct` = `Q_relax` / 블록 용량), 상대가 없으면 상대값을
   그대로 반환하는 fallback (`chgCapa_CCratio_norm`).
2. **실측 중복** — $\min(|r_{\text{pearson}}|, |r_{\text{spearman}}|)$의
   **per-cell 중앙값이 3셀 이상에서 $\ge 0.98$**. DOE1·DOE2·DOE3·ASSB 6셀 911행.
   pooled 상관은 교차확인용으로만 기록하고 **단독 근거로 쓰지 않는다** — pooling은
   셀 간 분산을 상관계수에 섞는다. `chg_temp_avg` ~ `dchg_temp_avg`가 pooled
   0.996인 것은 셀 간 온도차 때문이고 표본이 충분한 유일한 셀에서는 0.75다.
   반대로 pooled가 낮아도 병합을 막지 않는다 — `EoD_chgR_10s` ~ `_T25`는 같은
   이유로 pooled 0.54지만 셀 내부에서 온도보정은 단조변환이라 per-cell 1.0000이다.
3. **같은 측정량의 같은 구간** — 상관만으로는 부족하다. 의도적으로 만든 분해축
   (SOC 밴드, SOC 지점, 펄스 지속시간)에서 서로 다른 위치에 있는 두 컬럼은
   상관이 높아도 병합하지 않는다. 그 축은 열화를 **국소화**하기 위해 존재하므로
   병합은 중복 제거가 아니라 축 삭제다. 병합은 **한 위치 안에서만** 한다 — 같은
   밴드의 loop area와 peak \|dV\|, 같은 펄스의 원본과 온도보정값.

#### 폐기된 정확한 별칭 (8개)

컬럼 자체를 없앴다. 옛 테이블은 `apply_aliases`가 이름을 매핑한다. 6셀 911행
pre-change 추출에서 모든 항등식을 재구성 검증했다 (최대 상대편차 5.7e-16).

| 폐기 | 대체 | 항등식 | 검증 |
|---|---|---|---|
| `CI` | `CE` | $CI = 100 - CE$ | n=878, max\|diff\|=0 |
| `ocv_spread_compression` | `delta_ocv_spread_20_80` | 같은 값을 두 번 대입 | n=75, max\|diff\|=0 |
| `ocv_spread_slope` | `delta_ocv_spread_20_80` | 등간격 3점 fit → $-\Delta\text{spread}/60$ | n=75, 2.4e-19 |
| `ocv_block_id` | `block_id` | broadcast 복제 | n=75, max\|diff\|=0 |
| `ocv_V_inf_soc{80,50,20}` | `V_inf_rest_soc{80,50,20}` | 동일 `self_discharge_for_cycle` fit을 두 번 | n=25×3, max\|diff\|=0 |
| `R_SOC_diff_20_80` | `R_SOC_slope` | 등간격 3점 fit → $-60 \times \text{slope}$ | n=72, 1.8e-15 |

`CI` 폐기 시 `CI_per_hour`를 `CE`에서 직접 계산하도록 바꿨다. 시간 정규화된
형태는 `CE`와 pooled 0.90 / per-cell 0.96 으로 **독립이므로 유지해야 한다.**

#### 검토 후 기각한 병합

측정했으나 기준 미달이라 병합하지 않았다. 기록해 두지 않으면 다음 사람이 같은
병합을 다시 제안한다.

| 후보 | 실측 | 기각 이유 |
|---|---|---|
| `V_inf_rest_soc*` ~ `V_inf_est_soc*` | pooled 0.972 / 0.991 / 0.907 | 장시간 rest fit vs 펄스 회복 외삽 — 서로 다른 추정기. 3개 중 2개 미달 |
| `EoD_restV_60s` ~ `EoD_restV_end` | per-cell 0.976 | 방전측 완화가 충전측보다 느리다. 충전측 `EoC_restV_60s`는 0.998로 병합됨 |
| `chg_temp_avg` ~ `dchg_temp_avg` | pooled 0.996, per-cell 0.75 | pooled 값이 셀 간 온도차 아티팩트 |
| `CE` ~ `CI_per_hour` | pooled 0.90 | 사이클 시간 정규화가 새 정보를 만든다 |
| `tau_CV` ~ `Q_CV_at_Tref` ~ `I_inf_norm` | 데이터 없음 | 같은 CV fit의 서로 다른 파라미터. 픽스처에서 한 번도 채워지지 않아 중복 근거가 0 |
| `dchg_dVdQ_peak{1,2,3}_Q` ~ `dchgCapa` | pooled 0.27 / 0.55, per-cell 0.39 / 0.76 | `regime_singularity`의 옛 제외 목록이 "capacity proxy"라 단정했으나 **근거 없음** — indicator로 복원 |
| `hyst_area` ~ `hyst_max_dV_mid` 등 밴드 교차쌍 | per-cell 0.982 ~ 0.995 (6셀) | **기준 2는 통과하나 기준 3에서 기각.** global 적분은 mid 밴드를 구조적으로 포함하고, 교차쌍끼리 서로 모순된다 (global+mid와 high+mid를 병합하면 전 밴드가 한 컬럼으로 붕괴). §5.11 SOC 분해 히스테리시스의 존재 이유인 국소화 축이 사라지므로 4개 패밀리를 유지하고 잔여 중복을 감수한다. low 밴드는 \|r\| = 0.004로 애초에 독립 |

#### role 재분류 (자동 pool 오염 제거)

| 컬럼 | role | 근거 |
|---|---|---|
| `dchg_dVdQ_SOC0_Q` | `target` | `dchgCapa`와 min-\|r\| **1.0000** — 방전 종료점의 Q, 즉 용량 그 자체 |
| `dSoHQ_dN`, `d2SoHQ` | `target` | 용량 상관은 0.34 / 0.10로 낮지만 **타깃을 미분한 값**. anomaly 증거로 쓰면 순환논리 |
| `dchg_cliff_sg_width_ah` | `qc` | $q_{span} \times sg\_window / (n_{interp}-1)$ — Savitzky-Golay 창 폭. 전 셀에서 용량의 4.19 %로 고정 비례. 셀이 만든 값이 아니라 **처리 파라미터**이므로 capacity 패밀리도 아니다 |

반대로 registry는 구 prefix 필터가 **이름만 보고 버렸던** DC-IR 분해 계열
(`R_ohmic_*`, `R_ct_*`, `R_diff_frac_*`, `R_recovery_tau*`, `R_SOC_slope` 등 55개)을
pool에 새로 넣는다. 파이프라인이 내는 가장 직접적인 임피던스 증거인데
`EoC_`/`EoD_`/`chg_`/`dchg_`/`delta_`/`f_` prefix에 걸리지 않아 anomaly score가
보지 못하고 있었다.

#### mode weight 중복 증거 제거

한 mode의 evidence 목록에 같은 패밀리가 두 번 들어가면 그 물리량이 mode 점수를
두 번 지배한다. **패밀리 총 가중치는 보존**하고(fold) 대표 항만 남겼다.

| config · mode | 제거 | 유지 | 가중치 | 911행 점수 변화 |
|---|---|---|---|---|
| `assb_si_v1` · LLI | `CI` | `CE` | 1.0+0.9 → 1.9 | **0** (`decrease_from_100` 스케일이 동일해 항등) |
| `fullcell_v1` · LAM_NE | `chgCapa_CCratio` | `chgCVtime` | 0.6+0.7 → 1.3 | median +0.000, max +0.104 |
| `fullcell_v1` · impedance | `delta_hyst_max_dV` | `delta_hyst_area` | 0.8+0.7 → 1.5 | median −0.044, max −0.078 |

두 fullcell mode 모두 총 가중치가 보존된다 (6.7 → 6.7, 8.0 → 8.0). 점수가
움직이는 건 남은 대표항의 evidence 값이 제거된 항과 다르기 때문이다.

LAM_NE에서 `chgCVtime`을 남긴 이유: `chgCapa_CCratio`는 92–100 범위의 **절대
퍼센트**인데 `direction: decrease`, `scale: 5`로 채점돼 $\tanh(-CCratio/5)$가
전체 911행에서 −1.0이었다. positive evidence를 낼 수 없고 **영구적으로 자기
자신을 conflicting feature로 보고**하던 항이다 (별개의 기존 config 버그).

`evidence_count`는 mode당 1 줄어든다 (LAM_NE 8→7, impedance 10→9, LLI 4→3).
같은 측정량 두 개가 독립된 두 증거로 계상되지 않게 된 결과이므로 의도된 변화다.
LLI의 경우 `LLI_pattern_score`는 911행에서 **비트 단위로 동일**하고
`LLI_confidence`만 최대 0.075 움직인다.

#### 폐기한 ad-hoc 제외 목록

registry가 role·family를 대신 판정하므로 손으로 관리하던 목록 4개를 없앴다.

| 위치 | 대체 |
|---|---|
| `models/predict.py` `NUMERIC_PREFIX` + `_META_COLS` + prefix 화이트리스트 | `primary_indicator_columns` |
| `analysis/indicator_screen.py` `_META` | role 태깅 + `one_per_family` |
| `analysis/indicator_screen_plots.py` `_EXCLUDE_FROM_RANK` | role + `is_family_primary` |
| `analysis/regime_singularity.py` `_EXCLUDE` (18개) | `primary_indicator_columns` |
| `tools/run_lges_lifetime_ranking.py` `EXCLUDE` (14개) | registry + 누출 방지 `TARGET_DERIVED`(5개)만 잔존 |

#### 회귀 검증

동일 6셀·동일 사이클 선택으로 재추출: **911행 유지**, 393 → 385 컬럼
(= 폐기 별칭 8개, 그 외 소실 0, 신규 0). 잔존 385개 컬럼 중 값이 바뀐 것은
`LLI_confidence` / `LLI_evidence_count` / `diagnosis_quality_score` 3개뿐이고
모두 위 CE/CI fold의 의도된 결과다. `relax_completeness_soc*`는 채워진 행이
75 → 25로 줄지만 이는 블록 3사이클에 같은 값을 복제하던 broadcast가 사라진 것으로,
값이 다른 행은 없고 블록 단위 QC는 `relax_completeness_max`가 그대로 담당한다.

### 4.5 두 트랙 분리 — 지표 점수 vs 원인 해석

물리화학적 해석(LLI / LAM_NE / Si hysteresis / contact …)은 **맞다는 전제**를
유지하되, 실행 순서를 둘로 나눈다.

| 트랙 | 질문 | 모듈 | 기본 입력 |
|---|---|---|---|
| **A. 지표 점수** (선행) | 어떤 지표가 얼마나 움직였나? | `models/indicator_scoring.py` | **routine only** (RPT/post-RPT/DC-IR 제외), 패밀리당 대표 1개 |
| **B. 원인 해석** (후행, 분리) | 왜 움직였나? 어느 모드인가? | `diagnosis/*` + mode_weights | 트랙 A에서 신뢰된 증거만 사용 |

트랙 A 출력 (`score_layer=indicator`):
- `indicator_score` / `indicator_flag` / `indicator_top` (사이클 롤업)
- `indicator_summary` (지표별 점수, mode 라벨 없음)
- CLI: `python -m cyclediag.tools.export_indicator_scores …`

트랙 B 출력 (`*_pattern_score`, `degradation_mode` …)은 extract의
`with_diagnosis=True`로만 켜며, 지표 점수 테이블에 mode 이름을 섞지 않는다.

API:
- `api.score_dataframe` → 트랙 A만
- `api.diagnose_dataframe` → 트랙 A (+ 기술 스크린). 이름과 달리 **원인 진단이 아님**
- 원인 진단 → `cyclediag.diagnosis.diagnose_feature_table`

이후 개선(CE 강등, cliff 3분량, hyst 밴드 점수, Q_relax 정의 단일화 등)은
**트랙 A 점수 품질**과 **트랙 B 매핑 품질**을 따로 올린다.

### 4.6 두 트랙 후속 TODO — **구현됨**

물리화학 해석을 전제로 제안했던 P0–P3를 적용했다. 지표 점수(A)와 원인
진단(B)은 계속 분리한다.

| 우선순위 | 항목 | 적용 위치 |
|---|---|---|
| P0 | 진단도 routine-only 게이트 | `diagnosis/engine.py` → `filter_scoring_rows` |
| P0 | CE 강등, LLI에 OCV/VE/EoD 강화 | `mode_weights_*_v1.json` |
| P1 | cliff 3분량 + absolute-Ah | LAM_NE / contact_loss 가중치 |
| P1 | hyst low-band (`hyst_frac_low` / `hyst_area_low`) | contact / contact_loss; impedance에서 total hyst 제거 |
| P1 | Health / Mechanism / Anchor 리포트 층 | `analysis/indicator_layers.py`, `api.score_dataframe["by_layer"]` |
| P2 | `Q_relax_rpt` / `Q_relax_dcir` + primary + `Q_relax_source` | `features/enrich_assb.py` |
| P2 | R 성분·랜드마크 앵커 요약 | `analysis/resistance_anchors.py` |
| P2 | LLI vs kinetic-termination 판별 | `analysis/lli_kinetic_split.py` |
| P2 | SoHQ 구간 셀 비교 + knee split | `analysis/sohq_interval_compare.py` |
| P3 | temp / stack pressure / halfcell 신뢰도 캡 | `diagnosis/constraints.py` |
| P3 | set4 brief 제품 경로 | `tools/export_set4_brief.py` |

```bash
python -m cyclediag.tools.export_set4_brief \
  --fixtures example/fixtures/doe/DOE1/set4_SJ900 \
  --out example/output/set4_brief --with-causal
```

---

## 5. 알고리즘 명세

각 항목은 **목적 / 입력 / 알고리즘 / 출력 / 파라미터 / 실패조건 / 검증** 순으로 기술한다.
구현 모듈 경로를 함께 명시한다.

---

### 5.1 dQ/dV 필터 파라미터 스윕 진단

> **모듈:** `cyclediag/tools/diagnose_dqdv_filter_sweep.py`
> **우선순위: 1** — 이 결과가 §5.8(제약 디컨볼루션)의 필요 여부를 결정한다.

#### 목적

0.5C에서 dQ/dV 피크가 뭉치는 현상이

- **(A)** 실제 물리적 병합 (분극 이동 + broadening)
- **(B)** 스무딩/보간 파라미터에 의한 처리 인공물

중 무엇인지 판정한다. (B)라면 §5.8 전체가 불필요해진다.

#### 입력

- raw CSV
- `rpt_cycle` — C/3 RPT 사이클 (예: 107)
- `routine_cycle` — 0.5C routine 사이클 (예: 50)
- `leg` — `charge` | `discharge` (둘 다 실행 권장)

#### 알고리즘

```
1. 컬럼 해석
   논리명 ← 실제 컬럼명 매핑 (cycle/voltage/current/capacity)
   후보 리스트로 자동 탐색, 대소문자 무시

2. leg 추출 — step_type 라벨에 의존하지 않음
   cur ← Current
   active ← |cur| > rest_current_max        # 72Ah 셀이면 0.5 A 권장
   discharge_leg ← active & (cur < 0)
   charge_leg    ← active & (cur > 0)
   (전류 부호 규약이 반대인 장비 대비 fallback: 샘플 수 < 50이면 부호 반전 재시도)
   (전류 컬럼 부재 시: 용량 차분 dQ의 부호로 대체)

3. 파라미터 격자 순회
   for (n_interp, sg_window, sg_poly, merge_dv, min_distance_frac) in GRID:

     3a. dQ/dV 계산
         Q ← |Capacity - Capacity[0]|         # 누적 절대량으로 단조화
         (Q, V) ← Q 오름차순 정렬, 중복 Q 제거 (diff > 1e-9)
         Q_grid ← linspace(Q.min, Q.max, n_interp)
         V_grid ← interp(Q_grid, Q, V)
         win ← 홀수 보정, [poly+2, len//2*2-1]로 clip
         V_s ← savgol(V_grid, win, poly)      # smooth_then_diff
         Q_s ← savgol(Q_grid, win, poly)
         dQdV ← gradient(Q_s) / gradient(V_s)   # |dV|<1e-9는 NaN 후 보간
         dQdV ← |dQdV|

     3b. 물리 단위 환산 — 핵심 진단량
         ah_per_pt ← median(|diff(Q_grid)|)
         mv_per_pt ← median(|diff(V_s)|) * 1000
         smooth_Ah ← ah_per_pt * sg_window
         smooth_mV ← mv_per_pt * sg_window          ★ 이 값이 진단의 핵심
         min_dist_Ah ← ah_per_pt * min_distance_frac * n_interp

     3c. 피크 검출
         span ← max(dQdV) - min(dQdV)
         idx, props ← find_peaks(dQdV,
                          prominence = prominence_frac * span,
                          distance   = min_distance_frac * n_interp,
                          width      = min_width_points)
         # merge_dv 병합 (현행 로직 재현)
         for i in 1..len(idx)-1:
             if |V[idx[i]] - V[idx[last_kept]]| < merge_dv:
                 더 높은 쪽만 유지
             else: keep

     3d. 분해능 계산
         widths_pts ← peak_widths(dQdV, idx, rel_height=0.5)
         widths_V   ← widths_pts * median(|diff(V_grid)|)
         for 인접쌍 (k, k+1):
             Rs[k] ← |V[k+1]-V[k]| / (1.18 * (w[k] + w[k+1]))
         Rs_min, Rs_med ← min, median

     3e. 인공물 위험도
         smooth_vs_merge ← smooth_mV / (merge_dv * 1000)
         risk ← HIGH if ≥1.0 else MED if ≥0.4 else LOW

4. 판정
   rpt_ref ← median(rpt_n_peaks over valid combos)     # 전역 최대 아님
   recovered ← combos where (routine_n_peaks ≥ rpt_n_peaks at SAME params)
                        and (routine_n_peaks ≥ rpt_ref)
   if recovered 비어있지 않음 → 판정 (B) 처리 인공물
   else                       → 판정 (A) 실제 병합 → §5.8 필요

5. 권장 조합 정렬
   recovered.sort_by(Rs_min DESC, smooth_mV ASC)
   ※ Rs_min ≥ 1.5 인 조합 우선. Rs_min ~0.3 은 노이즈 과검출 의심
```

#### 파라미터 격자 (기본값)

```python
GRID = {
    "n_interp":          [500, 1500, 3000],
    "sg_window":         [7, 11, 21, 31],
    "sg_poly":           [3],
    "merge_dv":          [0.012, 0.006, 0.003],
    "min_distance_frac": [0.04, 0.015],
}
# 고정: prominence_frac=0.02, min_width_points=5
```

#### 출력

```
sweep_results.csv
  n_interp, sg_window, sg_poly, merge_dv, min_distance_frac,
  rpt_n_peaks, routine_n_peaks,
  {rpt,routine}_ah_per_pt, _mv_per_pt, _smooth_Ah, _smooth_mV,
  _min_dist_Ah, _Rs_min, _Rs_med,
  smooth_vs_merge, artifact_risk

sweep_summary.txt   현행 기준선 + 판정 + 권장 조합 상위 10
sweep_overview.png  (sg_window × n_interp) 피크개수 히트맵 + 곡선 오버레이
```

#### 실패 조건

| 조건 | 처리 |
|---|---|
| leg 샘플 < 50 | 예외. `--leg` 또는 `rest_current_max` 조정 안내 |
| 유효 (Q,V) < 20 | 해당 조합 스킵, `n_peaks = -1` |
| `span ≤ 0` | 피크 0개 반환 |
| 피크 < 2 | Rs = NaN (판정에서 제외) |

#### 검증

- 합성 곡선(알려진 피크 3개 + 인위적 broadening 2.2×)으로 스윕 → (B) 판정이 나오는지 확인
- 충전·방전 leg를 모두 실행 (Si 히스테리시스로 분해능이 다를 수 있음)
- `--list-cycles`로 전류 크기별 사이클을 먼저 확인 (0.33C ≈ 25.78 A, 0.5C ≈ 38.67 A, 1C 펄스 ≈ 77.34 A)

---

### 5.2 적응형 스무딩 (V-폭 기준)

> **모듈:** `cyclediag/features/dqdv_peaks.py` 수정
> **선행조건:** §5.1이 (B) 판정을 냈을 때 우선 적용. (A)여도 여전히 개선임.

#### 목적

Q축 등간격 보간 + 고정 인덱스 SG 필터의 조합은
**전압 영역별 실효 스무딩 폭을 제각각으로 만든다.**
평탄 구간에서는 과도하게 뭉개고, 급경사 구간에서는 거의 스무딩하지 않는다.

#### 알고리즘

```
1. 목표 스무딩 폭을 물리량으로 지정
   target_smooth_mV  (기본 4.0 mV — merge_dv 12 mV의 1/3)

2. 국소 전압 밀도 계산
   Q_grid ← linspace(Q.min, Q.max, n_interp)
   V_grid ← interp(Q_grid, Q, V)
   local_mv_per_pt[i] ← |V_grid[i+1] - V_grid[i-1]| / 2 * 1000

3. 국소 윈도우 크기
   win[i] ← clip( round(target_smooth_mV / local_mv_per_pt[i]) | 홀수화,
                  win_min = poly+2,
                  win_max = n_interp // 8 )

4. 가변 윈도우 SG 적용
   방법 A (권장, 단순): 윈도우를 몇 개 구간(bin)으로 양자화하고
                        각 bin에 대해 savgol을 전체 배열에 적용한 뒤
                        해당 bin 구간만 취해 이어붙임 (경계는 선형 blend)
   방법 B (정밀): 각 점마다 국소 다항 회귀 (LOESS 유사). 비용 큼

5. 미분
   dQdV ← gradient(Q_s) / gradient(V_s)

6. 진단 출력 (반드시 함께 기록)
   effective_smooth_mV_median / p05 / p95
   → 이 값이 전 구간에서 target 근처로 균일해졌는지 확인
```

#### 대안 — 더 단순한 처방

가변 윈도우 구현이 부담이면, **보간축을 V로 전환**하는 것만으로도 상당 부분 해결된다.

```
interp_axis = "V"
V_grid ← linspace(V.min, V.max, n_interp)   # 전압 등간격
Q_grid ← interp(V_grid, V_sorted, Q_sorted)
→ sg_window가 곧 일정한 mV 폭이 됨 (n_interp와 전압창으로 직접 환산)
   예: 전압창 1.7 V, n_interp=2000 → 0.85 mV/pt, sg_window=7 → 5.95 mV
```

단점: 급경사 구간(양 끝단)에서 Q 해상도가 떨어진다.
→ **Q축 보간은 dV/dQ 분석용, V축 보간은 dQ/dV 분석용**으로 분리하는 것이 정석이다.

#### 출력 컬럼 추가

```
dqdv_smooth_target_mV, dqdv_smooth_eff_mV_med,
dqdv_smooth_eff_mV_p05, dqdv_smooth_eff_mV_p95,
dqdv_interp_axis, dqdv_n_interp
```

#### 검증

- 합성 곡선에서 알려진 피크 폭이 보존되는지 (스무딩 전후 FWHM 비교)
- `effective_smooth_mV`의 p95/p05 비가 2 이내로 들어오는지
- 동일 사이클을 Q축·V축 양쪽으로 처리해 피크 개수·위치 일치도 비교

---

### 5.3 R(t) 3성분 분해

> **모듈:** `cyclediag/features/dcir_decompose.py` (신규)
> **우선순위: 2** — ASSB에서 접촉 손실 vs 계면 화학 열화를 가르는 유일한 저비용 수단.
> **전제 충족:** 10 Hz 전 구간 샘플링 + raw V–I 저장 확인됨.

#### 물리 모델

$$V(t) = V_0 - I\left[R_\Omega + R_{ct}\left(1 - e^{-t/\tau_{ct}}\right) + A\sqrt{t}\right]$$

| 성분 | 시간 영역 | ASSB 귀속 |
|---|---|---|
| $R_\Omega$ | $t \to 0$ | **접촉 면적 손실 · void** |
| $R_{ct}, \tau_{ct}$ | 1–10 s | 계면상 성장 · SE 분해 |
| $A$ (Warburg 유사) | 10–30 s | 고상 확산 (SE + 입자 내) |

#### 알고리즘

```
1. 펄스 구간 절단
   DC-IR 사이클(TC 4/5/6, 109/110/111, ...)에서
   pulse_start ← |I| 가 rest_current_max를 처음 초과하는 인덱스
   pulse_end   ← pulse_start + 30 s
   t ← time - time[pulse_start]      (0 ~ 30 s, 301점)
   V ← Voltage[pulse],  I ← median(Current[pulse])
   V0 ← Voltage[pulse_start - 1]     # 펄스 직전 마지막 rest 전압

2. 유효성 검사
   assert n_points ≥ 250            # 10 Hz × 30 s ≈ 301
   assert |I| > 0.5 * expected_pulse_current
   assert count(t ≤ 1.0) ≥ 8        # R_Ω 외삽에 필요
   assert std(I[1:]) / |median(I)| < 0.02   # 정전류 유지 확인

3. R(t) 곡선
   R(t) ← |V0 - V(t)| / |I| * 1000       [mΩ]

4. 단계적 피팅 (전역 최적화보다 안정적)

   4a. R_Ω — t→0 외삽
       초기 구간 t ∈ [0.0, 0.3] s 사용
       R(t) ≈ R_Ω + b·sqrt(t) 로 sqrt(t) 선형회귀 → 절편 = R_Ω
       (단순 R(0.1s) 채택보다 노이즈에 강함)
       fallback: 점이 부족하면 R(t_min) 사용 + flag

   4b. A (확산 계수) — 후반 구간
       t ∈ [10, 30] s 에서 R(t) vs sqrt(t) 선형회귀
       A ← 기울기,  R_late_intercept ← 절편
       (이 구간에서 지수항은 사실상 포화 → R_Ω + R_ct 로 수렴)

   4c. R_ct, tau_ct — 중간 구간 잔차
       resid(t) ← R(t) - R_Ω - A·sqrt(t)
       t ∈ [0.3, 10] s 에서
       resid(t) = R_ct (1 - exp(-t/tau_ct)) 를 비선형 최소자승
       초기값: R_ct0 = resid(10s), tau_ct0 = 2.0 s
       경계: R_ct ∈ [0, 5·R_Ω_scale], tau_ct ∈ [0.05, 20]

   4d. 전역 정련 (선택)
       위 4개를 초기값으로 전체 구간 curve_fit 1회
       조건수 > 1e8 이면 4c 결과 유지 (축퇴)

5. 품질 지표
   rmse ← RMS(R_fit - R_obs) / mean(R_obs)
   r2   ← 1 - SS_res/SS_tot
   cond ← 피팅 자코비안 조건수
   valid ← (rmse < 0.03) and (r2 > 0.98) and (cond < 1e8)

6. 파생 지표
   R_30s_total    ← R(30s)
   R_ohmic_frac   ← R_Ω / R_30s_total
   R_ct_frac      ← R_ct / R_30s_total
   R_diff_frac    ← A·sqrt(30) / R_30s_total
   (사이클 궤적에서)
   R_ohmic_growth_100 ← rolling slope of R_Ω per 100 cycles (% 기준)
   R_ct_growth_100    ← 동일
   mech_vs_chem_ratio ← R_ohmic_growth_100 / R_ct_growth_100
```

#### 출력 컬럼 (SOC별로 접미사 `_soc20/_soc50/_soc80`)

```
R_ohmic, R_ct, tau_ct, A_diff, R_30s_total,
R_ohmic_frac, R_ct_frac, R_diff_frac,
dcir_fit_rmse, dcir_fit_r2, dcir_fit_cond, dcir_fit_valid,
R_ohmic_growth_100, R_ct_growth_100, A_diff_growth_100,
mech_vs_chem_ratio
```

#### 해석 규칙 (진단 엔진에 등록)

| 조건 | 시사 모드 | 신뢰도 |
|---|---|---|
| `mech_vs_chem_ratio > 1.5` 지속 | **접촉 손실 지배** | 높음 (구속 압력 정보 있으면) |
| `R_ct_growth` ≫ `R_ohmic_growth` | 계면상 성장 / SE 분해 | 높음 |
| `A_diff_growth` 최대 | 고상 확산 제한 | 중간 |
| `R_ohmic_frac` 단조 상승 | void 형성 진행 | 중간 |
| `tau_ct` 증가 + `R_ct` 증가 | 계면층 **두께** 증가 | 중간 |
| `tau_ct` 감소 + `R_ct` 증가 | 유효 **면적** 감소 | 중간 |

마지막 두 줄이 유용하다. $R_{ct} \propto 1/(A \cdot i_0)$, $\tau_{ct} \approx R_{ct} C_{dl}$ 이고
$C_{dl} \propto A$ 이므로, 면적 감소는 $R_{ct}$↑ 와 $\tau_{ct}$ 거의 불변~감소를 낳고,
두께 증가는 둘 다 증가시킨다. **접촉 손실과 계면 성장의 구분 근거**가 된다.

#### 실패 조건

| 조건 | 처리 |
|---|---|
| t≤1 s 점 < 8 | `R_ohmic` NaN + `dcir_fit_valid=False` |
| 정전류 이탈 (std/median > 2 %) | 전체 무효 |
| `cond > 1e8` | 4d 스킵, 단계적 결과 유지 + flag |
| 1C 비선형 의심 | §5.3 주의 참조 — 절대값 대신 상대 추이만 사용 |

#### 1C 펄스 비선형성 주의

77.34 A ≈ 1C @ 45 °C는 과전압이 커서 Butler-Volmer 비선형 영역일 수 있다.

$$\eta_{ct} = \frac{RT}{\alpha F}\ln\frac{i}{i_0} \quad \text{(Tafel 영역)} \implies R_{ct} \text{ 가 전류 의존}$$

- **검증:** 0.5C 펄스를 1회 병행해 $R_{ct}$가 일치하는지 확인
- 불일치 시: $R_\Omega$는 여전히 유효(선형), $R_{ct}$는 **동일 조건 간 상대 추이로만** 사용
- 출력에 `dcir_linearity_verified: bool` 필드 추가

#### 검증

- 합성 트레이스(알려진 $R_\Omega, R_{ct}, \tau_{ct}, A$)로 복원 오차 측정 (목표 < 5 %)
- 동일 사이클 3 SOC 간 $R_\Omega$가 SOC 의존성이 작아야 함 (물리적으로 옴 성분은 SOC 무관)
  → 크게 변하면 피팅 또는 절단 구간 문제
- 장비 계산 `Impedance (ohm)` 컬럼과 $R(30s)$ 비교 (일치해야 정상)

---

### 5.4 자가방전율 / 미세단락 검출

> **모듈:** `cyclediag/features/self_discharge.py` (신규)
> **우선순위: 3** — 현재 완전히 미탐지인 실패 모드. **추가 실험 비용 0.**

#### 목적

ASSB의 주요 실패 모드인 **덴드라이트 관통 / 미세단락**을 검출한다.
DC-IR 각 펄스 전 rest 3600 s × 3 SOC가 이미 기록되어 있다.

#### 알고리즘

```
1. rest 구간 식별
   DC-IR 사이클에서 펄스 직전 rest 구간 (|I| ≤ rest_current_max, 지속 ≈ 3600 s)
   t ← 0 ~ 3600 s,  V(t)

2. 완화 성분과 분리 — 후반부만 사용
   확산 완화는 지수적으로 감쇠하므로 후반부에서는 무시 가능
   window ← t ∈ [1800, 3600] s
   선형회귀: V(t) = a + b·t
   self_discharge_rate ← -b * 3600 * 1000     [mV/h]

3. 완화 잔존 확인 (중요)
   완화가 아직 남아 있으면 자가방전으로 오인한다.
   전체 구간에 V(t) = V_inf - C·exp(-t/tau_relax) - k·t 를 피팅해
   tau_relax 를 추정하고,
   relax_residual_ratio ← C·exp(-1800/tau_relax) / |b·1800|
   if relax_residual_ratio > 0.3 → self_discharge 신뢰도 하향

4. 대안 추정 — 완화 모델 동시 피팅 (권장)
   V(t) = V_inf - C1·exp(-t/tau1) - C2·exp(-t/tau2) - k·t
   자유도 6. 3600 s에 10 Hz면 점이 충분하므로 식별 가능.
   self_discharge_rate ← k * 3600 * 1000  [mV/h]
   → 3의 잔존 완화 문제를 구조적으로 해결

5. 용량 환산 (선택)
   dQ/dV 를 이용해 mV/h → mAh/h 로 변환
   self_discharge_Q ← self_discharge_rate / 1000 * (dQ/dV @ 해당 SOC)
   [Ah/h] → 저항 환산: R_leak ≈ V / I_leak

6. SOC별 비교
   SOC 80에서 가장 민감 (전위차 최대)
   sd_ratio_80_20 ← rate(80) / rate(20)
```

#### 출력 컬럼

```
self_discharge_rate_soc20/50/80        [mV/h]
self_discharge_Q_soc20/50/80           [mAh/h]  (선택)
R_leak_est_soc80                       [ohm]    (선택)
sd_relax_tau1, sd_relax_tau2           [s]
sd_fit_valid, sd_relax_residual_ratio
sd_ratio_80_20
```

#### 해석 규칙

| 관측 | 시사 |
|---|---|
| `self_discharge_rate` 사이클에 따라 **단조 증가** | **관통 진행** — 최우선 경보 |
| SOC 80에서만 급증 | 고전위 구동 관통 / SE 산화 분해 |
| 전 SOC 균일 증가 | 화학적 자가방전 (shuttle) |
| 급격한 단발 상승 후 유지 | 부분 단락 발생 |

**경보 임계:** 초기 대비 3배 초과 또는 절대값 > 5 mV/h → `alert`

> 이 지표는 관통이 완전 단락으로 진행하기 전에 잡히는 **선행 지표**다.
> 용량 곡선에는 후기에나 나타난다.

#### 실패 조건

- rest 지속 < 1800 s → 산출 불가
- 온도 변동이 있으면 V 드리프트가 섞임 → **온도 로그 복구 후 재검토 필요** (현재 0.0)
- 3600 s에서도 완화가 지배적이면 (`tau_relax > 1200 s`) 신뢰도 하향

#### 검증

- 동일 셀의 3개 SOC rest에서 `tau_relax`가 유사해야 함 (물리적 일관성)
- BOL 셀의 자가방전율이 문헌 ASSB 값 범위(수 μA 급) 내인지
- 합성: 알려진 leak 전류를 주입한 rest 곡선에서 복원 확인

---

### 5.5 펄스 후 회복 τ

> **모듈:** `cyclediag/features/dcir_decompose.py` 내
> **전제 충족:** 펄스 후 rest 1800 s 기록 확인됨.

#### 목적

펄스 자체의 피팅(§5.3)과 **독립적인** 시정수 추정치를 얻어
확산 성분($A_{diff}$)을 교차검증한다. 피팅 축퇴를 잡아내는 안전장치.

#### 알고리즘

```
1. 회복 구간 절단
   펄스 종료 시점부터 rest 1800 s
   V_relax(t),  t = 0 ~ 1800 s
   V_end ← 펄스 마지막 전압

2. 2지수 피팅
   V(t) = V_inf - B1·exp(-t/tau_r1) - B2·exp(-t/tau_r2)
   초기값: tau_r1 = 5 s (전하이동 완화), tau_r2 = 300 s (확산 완화)
   경계: tau_r1 ∈ [0.5, 50], tau_r2 ∈ [50, 3000], tau_r1 < tau_r2

3. 진폭 비
   relax_amp_ratio ← B2 / (B1 + B2)     # 확산 성분 비중

4. 교차검증
   §5.3의 tau_ct 와 tau_r1 이 같은 자릿수여야 한다.
   |log10(tau_ct / tau_r1)| > 0.7  → 둘 중 하나가 축퇴 → 두 결과 모두 flag

   A_diff 와 B2 는 같은 방향으로 움직여야 한다.
   상관계수가 음수로 나오면 확산 항 해석 재검토
```

#### 출력 컬럼

```
R_recovery_tau1, R_recovery_tau2, relax_amp_ratio,
V_inf_est, recovery_fit_r2,
tau_consistency_flag           # tau_ct vs tau_r1 일치 여부
```

#### 부가 가치

`V_inf_est`는 해당 SOC의 **준-OCV 추정치**다.
사이클에 따른 `V_inf_est(SOC)` 이동은 반쪽셀 없이 얻는
**전극 슬리피지 / OCV 곡선 이동의 직접 증거**가 된다.

```
ocv_shift_soc50 ← V_inf_est_soc50(N) - V_inf_est_soc50(baseline)
```

SOC가 상대 % 기준이므로 이 값은 LLI와 LAM이 섞여 있으나,
§5.6의 곡선 정합과 결합하면 분리 가능하다.
---

### 5.6 3-파라미터 곡선 정합

> **모듈:** `cyclediag/features/curve_fit.py` (신규)
> **우선순위: 5** — 반쪽셀 없이 가능한 LAM/LLI/R 분해의 1차 근사.

#### 목적

`dchg_shape_DTW`는 "달라졌다"만 말한다.
곡선 변화를 **3개 물리 파라미터 + 잔차**로 분해한다.

$$V_N(Q) \approx V_{\text{ref}}\!\left(s \cdot Q + o\right) - I \cdot \Delta R$$

| 파라미터 | 물리 의미 | 모드 |
|---|---|---|
| `s` | 활물질 용량 축소 | LAM |
| `o` | 전극 슬리피지 / Li 재고 손실 | LLI |
| `ΔR` | 분극 증가 | 임피던스 |
| `residual` | 미설명 형상 왜곡 | 국부 열화, 불균일성 |

#### 알고리즘

```
1. 기준 곡선 준비
   baseline_cycle ← formation 이후 첫 RPT (현행 cycle=1 아님, §9 참조)
   ref: (Q_ref, V_ref) — 동일 leg, 동일 rate
   ※ Si 히스테리시스 때문에 충전은 충전 baseline, 방전은 방전 baseline

2. 대상 곡선
   (Q_N, V_N) — 동일 leg, 동일 rate만 비교 (0.5C↔0.5C, C/3↔C/3)

3. 목적함수
   V_model(Q) = interp(s·Q + o, Q_ref, V_ref) - I·dR
   loss(s,o,dR) = Σ_k w_k · [V_N(Q_k) - V_model(Q_k)]²

   가중 w_k: 양 끝단(급경사)은 낮게, 평탄부는 높게
   w_k ← 1 / (1 + (dV/dQ|_k / median(dV/dQ))²)
   → 끝단 노이즈가 피팅을 지배하는 것을 방지

4. 최적화
   초기값: s=1, o=0, dR=0
   경계:   s ∈ [0.5, 1.2], o ∈ [-0.3·Q_max, 0.3·Q_max], dR ∈ [-5, 50] mΩ
   방법:   Trust Region Reflective (scipy least_squares, loss='soft_l1')
           soft_l1 로 국부 이상치(피크 왜곡)의 영향 완화

5. 잔차 분석
   resid(Q) ← V_N(Q) - V_model(Q)
   fit_residual_rms  ← RMS(resid) * 1000        [mV]
   fit_residual_max  ← max|resid| * 1000        [mV]
   fit_residual_argmax_SOC ← Q(argmax|resid|) / Q_max * 100
   fit_r2

   ★ residual의 SOC 위치가 진단적이다:
     저 SOC 집중 → 음극 측 국부 열화
     고 SOC 집중 → 양극 측
     전 구간 분산 → 불균일성 / 다상 열화

6. 모드 귀속
   LAM_curve_proxy  ← (1 - s) * 100          [%]
   LLI_curve_proxy  ← o / Q_ref_max * 100    [%]
   R_curve_proxy    ← dR                     [mΩ]

   ⚠️ 이들은 proxy다. 전-셀 곡선만으로는 LAM_PE/LAM_NE 분리 불가.
      Layer 2 상태 지표로 등록하되 `*_est` 이름은 쓰지 않는다.

7. 축퇴 진단
   피팅 자코비안 J 에서 파라미터 상관행렬 계산
   corr(s, o) 가 |0.9| 초과면 두 파라미터 분리 불가 → 두 값 모두 flag
   → 이것이 "mode collision 감점"의 물리적 대체물이다
```

#### 출력 컬럼

```
{chg,dchg}_fit_scale, _fit_offset, _fit_dR,
_fit_residual_rms, _fit_residual_max, _fit_residual_argmax_SOC, _fit_r2,
_fit_corr_s_o, _fit_degenerate_flag,
LAM_curve_proxy, LLI_curve_proxy, R_curve_proxy
```

#### 검증

- 합성: ref 곡선에 알려진 (s, o, dR)을 적용 → 복원 오차 < 2 %
- `dR`이 §5.3의 `R_30s_total` 추이와 같은 방향인지 상관 확인
- `LAM_curve_proxy`가 §5.8 피크 면적 감소와 일관되는지

---

### 5.7 ΔQ(V) 통계

> **모듈:** `cyclediag/features/dqv_stats.py` (신규)
> **우선순위: 6** — 구현 난이도 대비 조기 예측력 최고.
> **참고 구현:** [microsoft/BatteryML](https://github.com/microsoft/BatteryML) Severson feature
> (초기·후기 방전곡선을 공통 전압축 보간 → ΔQ(V)의 var / skew / kurtosis / min).
> VP 진단 → **수명(RUL) 예측** 확장 시 1순위 벤치마크. 상세는 [§12.1](#121-batteryml에서-가져올-것).

#### 목적

Severson et al. (2019)의 핵심 feature. **용량 자체보다 훨씬 먼저 움직인다.**
LFP 124셀에서 초기 100사이클만으로 수명 예측을 가능케 한 지표.

#### 알고리즘

```
1. 고정 전압 그리드
   V_grid ← linspace(V_lo, V_hi, 1000)
   V_lo, V_hi ← 모든 비교 대상 사이클에서 공통으로 존재하는 전압 범위
                (양 끝 1 % 절사)

2. Q(V) 재표현
   각 사이클 방전 leg에서
   Q_N(V) ← interp(V_grid, V_sorted_desc, Q_sorted)
   ※ 방전은 V 단조 감소 → 정렬 방향 주의
   ※ 동일 rate끼리만 비교 (0.5C↔0.5C)
   ※ BatteryML과 동일: 사이클 간 공통 V축으로 보간 후 차분

3. 차이 곡선
   dQ_N(V) ← Q_N(V) - Q_ref(V)
   ref ← baseline_cycle (formation 후 첫 RPT 또는 지정 사이클)
   조기 예측 변형: Q_late(V) - Q_early(V)  (예: cycle 100 − cycle 10)

4. 통계량
   dQV_min      ← min(dQ_N)
   dQV_mean     ← mean(dQ_N)
   dQV_var      ← var(dQ_N)
   dQV_log_var  ← log10(var + eps)      ← Severson의 핵심 예측자
   dQV_skew     ← skew(dQ_N)
   dQV_kurtosis ← kurtosis(dQ_N)
   dQV_argmin_V ← V_grid[argmin(dQ_N)]  ← 손실 집중 전압 (진단적)

5. 짧은 구간 변형 (조기 예측용)
   Severson 원본은 cycle 100 - cycle 10 을 사용.
   본 프로토콜에서는 RPT 블록 간 차이도 유용:
   dQV_rpt_delta ← Q(V) at RPT_k - Q(V) at RPT_{k-1} 의 통계
```

#### 출력 컬럼

```
dQV_min, dQV_mean, dQV_var, dQV_log_var, dQV_skew, dQV_kurtosis,
dQV_argmin_V, dQV_ref_cycle, dQV_valid_V_range
dQV_rpt_delta_log_var        (RPT 간 변형)
```

#### 해석

- `dQV_log_var`의 초기 사이클 값이 **RUL과 강한 음의 상관**
- `dQV_argmin_V`가 이동하면 손실이 발생하는 전기화학 영역이 바뀐 것
- `dQV_skew` 부호는 손실이 고전압/저전압 어느 쪽에 치우쳤는지

#### 검증

- Ch22 (564 cycles, knee 포함 가능)에서 초기 100사이클 `dQV_log_var` vs 실제 수명 상관
- 공개 데이터셋 Severson으로 재현 확인 (원 논문 값과 비교 가능)

---

### 5.8 Q-domain 제약 디컨볼루션 + 피크 그룹

> **모듈:** `cyclediag/features/peak_deconv.py` (신규), `peak_assign.py` 확장
> **선행조건: §5.1이 (A) 판정을 냈을 때만 구현한다.**
> (B)면 §5.2 필터 수정으로 충분하므로 이 절 전체가 불필요하다.

#### 왜 V-domain 매칭이 실패하는가

0.5C에서의 병합은 **강체 평행이동이 아니다.**

| 원인 | 효과 | 균일성 |
|---|---|---|
| 옴 분극 $IR_\Omega$ | 균일한 V 이동 | 균일 |
| 전하이동 과전압 | 반응별 $i_0$ 에 따라 상이 | 비균일 |
| 고상 확산 제한 | 폭 확장 + **비대칭 tailing** | 비균일 |
| 반응 불균일성 | broadening | 비균일 |

따라서 ΔV만 빼는 보정은 원리적으로 실패한다.

#### 핵심 전략 — Q축에서 어사인

상전이 순서는 열역학적으로 결정되므로 전이는 **거의 동일한 통과 전하량**에서 일어난다.
분극은 V축을 압축·이동시키지만 Q축은 (양 끝 kinetic 손실 제외) 상대적으로 보존된다.

```
1차 공간: dV/dQ vs Q_norm     ← 어사인은 여기서
2차 공간: dQ/dV vs V          ← 확인용
```

#### 알고리즘

```
1. 앵커 확보 (RPT)
   for each RPT block (TC 4-6, 109-111, 214-216, ...):
       ★ cycle 2를 앵커로 사용 (cycle 1은 직전 0.5C 이력 보유)
       peaks_rpt ← detect_peaks(RPT_cycle2, domain="Q_norm")
       Q_norm ← Q / Q_total    (kinetic 손실 보정은 §5.6의 s, o 활용)
       anchor[block] ← {K, Q_k, A_k, sigma_k}

   noise_floor ← |RPT_cyc1 - RPT_cyc2| 피크 위치·높이 차이
                 (확정값: 용량 기준 0.065 %)

2. 앵커 간 보간 (RPT 주기 105 사이클 → 보간 구간이 길다)
   for routine cycle n between anchor[i] and anchor[i+1]:
       prior_Q_k(n) ← 선형 또는 spline 보간
       prior_A_k(n) ← 동일
   ※ 105 사이클 보간은 리스크가 크다 → 4의 검증 필수

3. 제약 디컨볼루션 (핵심)

   모델:
     (dQ/dV)_0.5C(V) ≈ Σ_k A_k · g(V; V_k^prior + Δ0 + β·Q_k, sigma_k, skew)
     g = EMG (exponentially modified Gaussian)  ← 확산 tailing 표현

   자유도 관리 (성패를 가름):
     K          : RPT에서 고정 (자유 파라미터 아님)
     A_k        : prior 대비 ±tol (기본 30 %) 내로 제약
     위치       : 전역 Δ0 + Q 선형항 β  ← 피크별 자유 이동 금지
     폭         : sigma_k = sqrt(sigma_k_prior² + sigma_rate²)
                  sigma_rate는 전역 1개  ← 피크별 자유 금지
     비대칭     : 전역 skew 1개

     → 자유 파라미터 3~5개
     ✗ 8피크 × 3파라미터 = 24개 자유 피팅은 절대 식별 불가

   최적화:
     least_squares(loss='soft_l1'), 초기값 = prior
     정규화항: λ·Σ(A_k - A_k_prior)² 추가

4. 순차 베이지안 (시간 연속성)
   사이클 n의 posterior를 사이클 n+1의 prior로 전달
   → 병합 구간의 모호성이 시간 연속성으로 해소된다
   Kalman smoother 또는 단순 EMA로 구현 가능

5. 분해능 판정 → 그룹화
   for 인접쌍 (k, k+1):
       Rs ← |V_{k+1} - V_k| / (1.18·(w_k + w_{k+1}))
       if Rs < 0.8:  merge_into_group(k, k+1)

   그룹 단위 물리량만 보고 (개별 분리 불가여도 강건):
       group_area     ← Σ A_k                 [매우 강건]
       group_centroid ← Σ A_k·V_k / Σ A_k     [강건]  (Q 버전도)
       group_width    ← 2차 모멘트            [강건]
       group_skew     ← 3차 모멘트            [중간]
   개별 peak_k_V/H는 NaN + 사유 코드로 보고

6. 식별성 진단
   cond ← condition_number(Jacobian)
   param_corr ← 파라미터 상관행렬
   peak_confidence ← f(cond, residual_rms / noise_floor, Rs)
```

#### 출력 스키마

```
peak_assign_method    : "rpt_anchored_deconv" | "direct_detect" | "group_only"
peak_k_resolved       : bool
peak_k_group_id       : int
peak_k_confidence     : 0–1
peak_k_V, peak_k_Q, peak_k_H, peak_k_sigma   (resolved=True 일 때만)
group_g_area, group_g_centroid_V, group_g_centroid_Q,
group_g_width, group_g_skew
deconv_residual_rms, deconv_cond, deconv_n_free_params
sigma_rate_global, shift_global, skew_global
```

#### 파생 지표

```
peak_k_sigma_growth_100     # 반응 불균일성, 입자 균열
group_g_width_growth_100
peak_resolution_Rs          # C/3에서의 Rs 하락 = 순수 재료 열화 신호
Q_PE_apparent               # 양극 특성 피크 간 Q 간격 (LAM_PE)
                            # ※ Si 음극이므로 NE 버전은 산출하지 않음
```

#### 검증 (반드시 구현)

| 방법 | 내용 |
|---|---|
| **Leave-one-RPT-out** | anchor $n_1$ + 0.5C 궤적으로 $n_2$ RPT 피크 예측 → 실측 비교. **유일한 준-정답 검증** |
| **RPT 2사이클 반복성** | cycle 1 vs 2 차이 = 노이즈 하한. 이보다 작은 이동은 무의미 |
| **합성 데이터** | 알려진 피크 → 인위적 broadening/shift → 복원율 측정 |
| **잔차 검사** | residual > noise_floor 이면 K 또는 모델 형태가 틀림 |
| **조건수** | > 1e8 이면 축퇴 → 그룹 모드로 강등 |

**RPT 주기 105 사이클 리스크:** leave-one-RPT-out 오차가 크면
중간에 경량 RPT(C/3 1사이클)를 추가하는 프로토콜 변경을 검토한다.

---

### 5.9 η(SOC) × DC-IR 결합 R_eff

> **모듈:** `cyclediag/features/overpotential.py` (신규)
> **우선순위: 7**

#### 목적

두 rate 곡선(C/3 RPT, 0.5C routine)에서 **EIS 없이 SOC 전 구간 분극 프로파일**을 얻는다.
DC-IR 3점만으로는 얇지만, η(SOC)의 연속 형상과 결합하면 전 구간 보정이 가능해진다.

#### 알고리즘

```
1. 과전압 곡선
   동일 Q(또는 Q_norm)에서 두 곡선의 전압차
   eta(Q) ← V_C3(Q) - V_0.5C(Q)          [V]
   ※ 방전 기준. 충전은 부호 반대

   대상 사이클 선택: RPT 블록에 가장 가까운 routine 사이클
                     (POST_RPT_EXCLUDE=5 를 지나서)

2. SOC 격자로 리샘플
   SOC ← (1 - Q/Q_total) * 100   (방전 기준)
   eta_grid ← interp(SOC_grid=[0..100], SOC, eta)

3. 대표값
   eta_SOC20, eta_SOC50, eta_SOC80
   eta_max, eta_argmax_SOC        ← 제한 전극 식별의 핵심
   eta_mean, eta_slope_lowSOC (SOC 0-30 기울기)

4. DC-IR 3점과 결합 → 전 SOC R_eff
   shape(SOC) ← eta_grid(SOC) / eta_grid(50)         # 형상, 무차원
   scale      ← argmin_c Σ_{s∈{20,50,80}} [c·shape(s) - R_DCIR(s)]²
                (최소자승 1파라미터)
   R_eff(SOC) ← scale · shape(SOC)

   ★ 이것이 §5.8 pseudo-OCV 보정의 입력이 된다

5. 잔차 = 정보
   resid_s ← scale·shape(s) - R_DCIR(s)   for s in {20,50,80}
   잔차가 크면 30 s DC-IR과 0.5C 정상상태의 물리가 다르다는 신호
   → PER (§5.10)로 정량화
```

#### 출력 컬럼

```
eta_SOC20, eta_SOC50, eta_SOC80, eta_max, eta_argmax_SOC,
eta_slope_lowSOC, eta_mean,
Reff_scale, Reff_shape_fit_r2, Reff_resid_soc20/50/80
```

#### 해석 규칙

| 관측 | 시사 |
|---|---|
| **저 SOC에서 η 증가** | 음극 수송 제한 / Si 열화 |
| **고 SOC에서 η 증가** | 양극 측 제한 / PE 계면 |
| `eta_argmax_SOC` **이동** | **제한 전극의 전환** — 매우 진단적 |
| η 전 구간 균일 증가 | 옴 성분 지배 (접촉 손실) → §5.3 `R_ohmic`과 교차확인 |

#### DC-IR SOC 순서 효과 보정

3 SOC가 **3개 연속 사이클에 분산**되어 있으므로 후속 SOC는 앞선 펄스 이력을 갖는다.

```
보정: 매 블록에서 순서가 동일한지 확인 (80 → 50 → 20)
      순서 효과를 상수 offset으로 흡수:
      R_corrected(s) = R_raw(s) - c_order[position(s)]
      c_order는 전 블록에서 공통 추정 (블록 간 상대 비교에는 상쇄됨)
```

#### 검증

- `Reff_shape_fit_r2` > 0.9 이면 결합이 유효
- `eta_argmax_SOC`가 §5.6 `fit_residual_argmax_SOC`와 일치하는지 (독립 경로 교차검증)

---

### 5.10 Q_relax · RCF · PER

> **모듈:** `cyclediag/features/rpt_metrics.py` (신규)
> **우선순위: 4 (Q_relax는 즉시 가능)**

#### 5.10.1 Q_relax — RPT 2사이클 차이

```
Q_relax ← Q_dchg(RPT_cycle2) - Q_dchg(RPT_cycle1)
Q_relax_pct ← Q_relax / Q_dchg(RPT_cycle2) * 100

노이즈 하한 (확정): 0.045 / 68.9 ≈ 0.065 %
→ |Q_relax_pct| > 0.065 % 인 경우만 유의로 판정
```

**해석:** 사이클에 따라 `Q_relax`가 **증가**하면
농도 구배 축적 / 불균일 SOC 심화 → ASSB에서는 **접촉 불균일 심화** 신호.
액체 셀의 wetting 저하 해석은 적용하지 않는다.

#### 5.10.2 rate capability fade (RCF)

```
RCF(N) ← Q_0.5C(N) / Q_C3(nearest RPT)
```

| 거동 | 해석 |
|---|---|
| RCF 하락 + Q_C/3 유지 | **순수 kinetic 열화** — 계면/접촉/수송. 활물질은 생존 |
| RCF 유지 + Q_C/3 하락 | **열역학적 용량 손실** — LAM/LLI |
| 둘 다 하락, RCF가 빠름 | 복합, kinetic 우세 |

**이 한 지표가 Layer 2 상태를 1차 분리한다.** 가중합보다 물리적 근거가 명확하다.

```
파생: RCF_slope_100  (100사이클당 변화율)
      RCF_at_SoHQ90  (SoHQ 90 % 시점의 RCF — 셀 간 비교용 정규화)
```

#### 5.10.3 polarization excess ratio (PER)

```
dV_observed ← V_C3(Q) - V_0.5C(Q) at SOC 50    (= eta_SOC50)
dI          ← |I_0.5C - I_C3|                   ≈ 38.67 - 25.78 = 12.89 A
PER ← dV_observed / (dI · R_DCIR_50 / 1000)
```

| PER | 해석 |
|---|---|
| ≈ 1 | 분극이 30 s DC-IR로 설명됨 (옴 + 전하이동 지배) |
| **> 1 (정상)** | 30 s가 못 잡는 장시간 확산 분극 존재 |
| **증가 추세** | 수송 열화가 계면 열화보다 빠름 → 고상 확산 제한 진행 |

**PER은 §5.9의 `Reff_scale` 보정계수 α와 같은 물리량이다.**
보정 파라미터가 그대로 진단 지표가 되는 구조.

#### 5.10.4 RPT 회복 용량 분해

```
Q_recovery ← Q_dchg(RPT_cyc2) - Q_dchg(last 0.5C before RPT)
           = rate 성분 + relaxation 성분
rate 성분        ≈ Q_C3 · (1 - RCF)
relaxation 성분  ≈ Q_relax (5.10.1)
→ 잔차가 있으면 미설명 회복 → 추가 조사 대상
```

#### 출력 컬럼

```
Q_relax, Q_relax_pct, Q_relax_significant,
RCF, RCF_slope_100, RCF_at_SoHQ90,
PER, PER_slope_100,
Q_recovery, Q_recovery_rate_part, Q_recovery_relax_part, Q_recovery_residual
```

#### 검증

- `Q_relax`의 블록 간 분포가 노이즈 하한 근처에 모이는지 (BOL에서)
- `PER`과 §5.3 `A_diff`가 같은 방향으로 움직이는지 (독립 경로 교차검증)
- `RCF`와 §5.9 `eta_mean`의 음의 상관 확인
---

### 5.11 SOC 분해 히스테리시스

> **모듈:** `cyclediag/features/lges_extra_indicators.py` 확장
> **우선순위: 8** — Si-rich 음극이므로 **Si 열화의 직접 지표**.

#### 목적

`hyst_area` 단일 값은 Si 열화의 SOC 국소성을 감춘다.
Si는 특정 리튬화 구간에서 우선적으로 열화하므로 구간 분해가 필수다.

#### 알고리즘

```
1. Q_norm 정렬
   충전: Q_chg_norm  ← Q_chg / Q_chg_total
   방전: Q_dchg_norm ← Q_dchg / Q_dchg_total
   ※ 용량이 다르므로 각자 정규화. 절대 Q로 맞추면 왜곡

2. 공통 격자
   x ← linspace(0.02, 0.98, 500)      # 양 끝 2 % 절사 (급경사 배제)
   V_chg(x), V_dchg(x) ← interp

3. 구간별 면적
   bands = {"low": (0.00, 0.20), "mid": (0.20, 0.80), "high": (0.80, 1.00)}
   for name, (a, b) in bands:
       hyst_area_{name} ← ∫_a^b [V_chg(x) - V_dchg(x)] dx      [V]
       hyst_max_dV_{name} ← max over band
       hyst_argmax_{name} ← x at max

4. 전체 및 비율
   hyst_area_total
   hyst_frac_low  ← hyst_area_low / hyst_area_total
   hyst_frac_high ← hyst_area_high / hyst_area_total

5. 궤적
   hyst_area_{band}_growth_100
   hyst_frac_low_shift ← hyst_frac_low(N) - hyst_frac_low(baseline)
```

#### 해석 규칙

| 관측 | 시사 |
|---|---|
| `hyst_area_low` 우선 증가 | **Si 리튬화 초기 구간 열화** (Si는 저 SOC에서 활성) |
| `hyst_area_high` 우선 증가 | 양극 측 또는 고전위 계면 |
| `hyst_frac_low` 상승 추세 | Si 기여도 증가 → 팽창/수축 스트레스 누적 |
| 전 구간 균일 증가 | 옴 성분 (접촉 손실) → §5.3 `R_ohmic`과 교차확인 |

**주의:** Si의 히스테리시스는 열화 없이도 크다. **절대값이 아니라 추세**로만 판정한다.

---

### 5.12 페이드 지수 · knee 검출

> **모듈:** `cyclediag/analysis/sohq_inflection.py` 확장
> **선행 모듈:** `cyclediag/analysis/rpt_recovery.py` — RPT 회복 bump 보정
> **진단 도구:** `cyclediag/tools/diagnose_rpt_recovery.py`
> **검증 데이터:** Ch22 (564 cycles, SoHQ ~65 %) — knee 포함 가능성 높음

#### 5.12.0 RPT 회복 bump 보정 (knee 선행 단계)

**본체 현상 (0.5C끼리 비교):** RPT를 시리즈에서 빼도 남습니다. 복귀한 0.5C routine
사이클의 용량이 RPT **직전** routine보다 높게 나옵니다 — rate(C/3 vs 0.5C) 차이가
아니라 **RPT 블록 동안의 긴 정지** 때문입니다. `POST_RPT_EXCLUDE=5`는 이 현상을
완화하려는 것이지만, 고정 5사이클은 후기 λ 증가 시 과소보정됩니다.

**가역 / 비가역 분해:**

$$Q_{\text{obs}}(n) = Q_{\text{irreversible}}(n) + Q_{\text{reversible}}(n)$$

RPT 정지가 $Q_{\text{reversible}}$을 0으로 리셋 → 복귀 직후 bump.
- **A** = 정지 직전 축적된 가역 손실 (열화 상태 직접 측정)
- **λ** = 가역 손실 재축적 속도
- knee 검출 입력 = **$Q_{\text{irreversible}}$** (`SoHQ_corrected`) — A·λ는 **기록 후** 제거

| 후보 메커니즘 | 되돌아가는 것 | ASSB |
|---|---|---|
| 농도 구배 완화 | Li 재분배 | 공통 |
| SOC 불균일 해소 | 국부 과·저리튬 균질화 | 공통 |
| 기계적 크리프 | Si 팽창 계면 부분 재접촉 | **전고체 특유** |
| 계면 재정렬 | SE\|활물질 접촉 회복 | **전고체 특유** |

BOL에는 가역 손실 축적이 없어 bump ≈ 0. **후기에만 나타나는 bump**는 접촉/불균일
심화의 지표.

**처리 순서:**

```
1. SoHQ_routine (0.5C only) — knee 1차 입력
2. pre-RPT anchor (오염 없는 기준선)
   anchor_k = mean(SoHQ_routine[block_start-W : block_start])   # W=5, routine만
   Ch22 → 6점 anchor 시리즈 (knee 기준선 + RPT 검증)
3. 회복 모델 (고정 제외 대신)
   Q(n) = Q_trend(n) + A·exp(-(n - n_block_end)/λ)
4. bump_contamination (knee 왜곡 정량화)
   fade_rate_intra  ← inter-RPT routine 구간 기울기 (bump 포함)
   fade_rate_inter  ← 연속 anchor 간 기울기
   bump_contamination = fade_rate_intra / fade_rate_inter - 1
5. 선행 지표
   bump_onset ← A > Q_relax 노이즈 하한 (0.065 %) 첫 블록
   knee_onset ← d²SoHQ_corrected/dN² > 3σ (NaN 갭 유지)
   가설: bump_onset_cycle < knee_onset_cycle
```

**회복 지표 (제거 대상이 아니라 ASSB 불균일성 측정치):**

| 컬럼 | 의미 | 단위 |
|---|---|---|
| `rpt_recovery_amplitude` | 회복 진폭 A (가역 손실) | % SoHQ |
| `rpt_recovery_decay_cycles` | 회복 시정수 λ | cycles |
| `pre_rpt_anchor_sohq` | 블록 직전 anchor | % SoHQ |
| `bump_contamination` | intra/inter fade 비 | — |
| `bump_onset_cycle` | A > 0.065 % 첫 블록 | cycle |
| `permanent_step_pct` | 비가역 잔차 (c) | % SoHQ |

§5.10 `Q_relax`와 동일 물리(농도 구배 완화)를 블록 시간 스케일에서 본다. 두 지표가
같은 방향이면 교차검증.

**진단 (Ch22, 6 RPT 블록):**

```bash
python -m cyclediag.tools.diagnose_rpt_recovery \
  --stepemd "example/docs/features/M01Ch022/*Ch22*stepend.csv" \
  --out example/output/rpt_recovery_ch22
```

- 모든 블록을 상대 사이클 `(n − n_block_end)`로 겹쳐 그림 → λ 적절성, λ(수명) 추세
- `SoHQ_mixed` vs `SoHQ_routine` vs `SoHQ_corrected` 분리 플롯
- d²SoHQ 추정 시 NaN 갭 유지 또는 추세 보간 + 신뢰도 가중 (갭을 0으로 메우지 않음)

#### 5.12.1 페이드 지수

$$Q(N) = 1 - a N^{b}$$

```
1. rolling window (기본 100 사이클, stride 20)
   window 내 (N, SoHQ/100) 사용
2. 로그 변환 선형회귀
   log(1 - Q) = log(a) + b·log(N)
   ※ Q > 0.999 인 초기 구간은 제외 (log 발산)
3. b 추정 + 표준오차
   fade_exponent_b, fade_exponent_b_se, fade_fit_r2
```

| b 값 | 지배 기구 (ASSB 맥락) |
|---|---|
| ≈ 0.5 | 확산 제한 계면상 성장 (LLI 지배) |
| ≈ 1.0 | LAM / 접촉 손실 선형 진행 |
| **> 1.5** | knee — 급격한 접촉 손실, SE 분해 가속, 관통 |

#### 5.12.2 knee 검출

**입력:** §5.12.0 보정 후 `SoHQ_corrected` (routine + 회복분 제거). RPT/C/3 용량은
`SoHQ_rpt`로만 검증에 사용.

두 방법을 병행하고 일치 여부를 신뢰도로 사용한다.

```
방법 A — Bacon-Watts 변화점 모델
   SoHQ(N) = a0 + a1(N - N_k) + a2(N - N_k)·tanh((N - N_k)/gamma)
   자유도 5 (a0,a1,a2,N_k,gamma)
   N_k ← knee_cycle
   gamma ← 전이 폭 (작을수록 급격)

방법 B — bisector 방식
   초기 구간 선형 적합 L1, 후기 구간 선형 적합 L2
   두 직선의 각이등분선과 SoHQ 곡선의 교점 → knee_cycle
   (Fermín-Cueto et al. 방식)

일치 판정:
   |knee_A - knee_B| / total_cycles < 0.1 → knee_confidence = high

knee onset (곡률 시작점):
   knee_onset_cycle ← d²SoHQ/dN² 가 노이즈 3σ를 처음 초과하는 지점
   (Savitzky-Golay 2차 미분, window는 사이클 수의 10 %)

severity:
   knee_severity ← |slope_after - slope_before| / |slope_before|
```

#### 출력 컬럼

```
# §5.12.0 회복 보정
SoHQ_routine, SoHQ_rpt, SoHQ_corrected, cycle_role,
rpt_recovery_amplitude, rpt_recovery_decay_cycles,
pre_rpt_anchor_sohq, bump_contamination, bump_onset_cycle,
permanent_step_pct,

# §5.12.1–2 knee / fade
fade_exponent_b, fade_exponent_b_se, fade_fit_r2,
knee_cycle_A, knee_cycle_B, knee_cycle, knee_confidence,
knee_onset_cycle, knee_severity, knee_gamma,
dSoHQ_dN, d2SoHQ    (기존, NaN-aware on RPT gaps)
```

#### 검증

- Ch22에서 두 방법의 knee가 일치하는지
- Severson 공개셋에서 knee 검출률·위치 오차 비교

---

### 5.14 신규 지표 계열 (가역/비가역 분해 이후)

> 이번 세션에서 도출. 기존 RCF·η(SOC)·R 3성분·자가방전·Q_relax와 **중복 제외**.
> **Anomaly 입력:** A·B 계열은 `Q_reversible_fraction` 하나만; 나머지는 진단·스크리닝 전용.

#### A. 가역/비가역 분해 (`rpt_recovery.py`)

| 지표 | 정의 | anomaly |
|---|---|---|
| `Q_reversible_fraction` | $A_k / (100 - \text{SoHQ}_k)$ | **유일 입력** |
| `Q_rev_accumulation_rate` | $A_k / \lambda_k$ | 스크리닝 |
| `rest_sensitivity` | SOC별 pre-RPT rest 회복 기여도 | 스크리닝 |

#### B. 블록 내 궤적 (`rpt_recovery.py` 확장)

| 지표 | 정의 |
|---|---|
| `block_curvature` | inter-RPT SoHQ 2차항 |
| `block_shape_consistency` | 상대인덱스 overlay 상관/DTW |
| `block_end_slope` | 회복 소진 후 기울기 (= 순수 비가역 fade) |

#### C. quasi-OCV drift (`ocv_drift.py`) — **구현됨**

DC-IR 3600 s pre-pulse rest → SOC 80/50/20 세 점 $V_\infty$.

| 컬럼 | 의미 |
|---|---|
| `V_inf_rest_soc{80,50,20}` | pre-pulse quasi-OCV. 해당 SOC 사이클 행에만 존재 |
| `ocv_spread_20_80` | $V_{20} - V_{80}$ |
| `delta_ocv_*` | BOL 블록 대비 drift. 블록 3사이클 전체에 broadcast |
| `ocv_parallel_shift` | 평행 이동 → **LLI + kinetic termination proxy** (R30 상관 필수) |
| `delta_ocv_spread_20_80` | spread 변화 → **electrode imbalance proxy** (비대칭 LAM/정렬) |
| `ocv_drift_mode` | `lli_parallel` / `lam_spread` / `local_soc*` |
| `relax_completeness_soc*` | 3600 s 후 잔여 완화 비율. 해당 SOC 사이클 행에만 존재 |
| `relax_completeness_max` | 블록 내 최악값. 블록 3사이클 전체에 broadcast |

모듈: `cyclediag/features/ocv_drift.py` · `enrich_assb` 자동 attach.

§4.4 통폐합으로 `ocv_spread_compression`(= `delta_ocv_spread_20_80`),
`ocv_spread_slope`(= $-\Delta\text{spread}/60$), `ocv_V_inf_soc*`
(= `V_inf_rest_soc*` 동일 fit), `ocv_block_id`(= `block_id`) 4개 컬럼이 폐기됐다.

#### D–F (미구현)

- `Q_irrev_heat` / `heat_per_Ah` — ΔE/Ah 발열 프록시
- `step_start_curvature` — leg 시작 √t 형상 (mini DC-IR)
- `dcir_order_effect` — 80→50→20 순서 이력

**우선순위:** 1 `Q_reversible_fraction` · 2 `ocv_drift_across_soc` · 3 `block_end_slope` · …

---

### 5.13 데이터 품질 게이팅

> **모듈:** `cyclediag/features/quality.py` (신규)
> **현행 `diagnosis_quality_score` 를 대체한다.**

#### 현행의 문제

```
mode data_quality = n_available_evidence / n_configured_terms
```
"증거 컬럼이 채워졌는가"만 본다. 실제 데이터 품질과 무관하다.

#### 신규 알고리즘

```
사이클마다 raw에서 직접 계산:

1. samples_per_mV
   dQ/dV 유효 구간에서 |dV| 총합 / 샘플 수의 역수
   → 피크 검출 신뢰도 결정
   기준: ≥ 0.5 samples/mV 권장

2. v_noise_sigma
   짧은 구간(50점) rolling에서 선형 추세 제거 후 잔차 std
   → 피크 높이 오차 전파의 근원
   부수 효과: ADC 양자화 폭 역추정 가능 (#17 UNKNOWN 해결)
   quant_step ← min(|diff(V)| where diff != 0)

3. dqdv_snr
   median(peak_height) / v_noise_sigma 기반 dQ/dV 노이즈
   기준: ≥ 10

4. rest_sufficiency
   rest_duration / tau_relax_est
   기준: ≥ 3 (tau 추정 유효 조건)

5. pulse_sample_present
   DC-IR 펄스 t ≤ 1 s 샘플 수 (현행 11 → 충족)
   pulse_current_stability ← std(I)/|median(I)|

6. leg_completeness
   충·방전 leg의 V 범위가 기대 창의 몇 %를 덮는지
   중단·resume으로 잘린 사이클 검출

7. temperature_available
   온도 컬럼이 실제 값인지 (현재 전부 0.0 → False)

집계:
   quality_score ← 가중 기하평균 (0 하나면 전체 0이 되도록)
   quality_score = Π (clip(q_i / q_i_target, 0, 1)) ^ w_i

게이팅:
   for each 지표 group:
       if 필요 품질 미충족 → 해당 지표 NaN + reason code
   예:
       peak_* 는 dqdv_snr ≥ 10 이고 samples_per_mV ≥ 0.5 일 때만 산출
       tau_* 는 rest_sufficiency ≥ 3 일 때만
       R_ohmic 은 pulse t≤1s 샘플 ≥ 8 일 때만
```

#### 출력 컬럼

```
samples_per_mV, v_noise_sigma, quant_step_est, dqdv_snr,
rest_sufficiency, pulse_sample_count_1s, pulse_current_stability,
leg_completeness, temperature_available,
quality_score, quality_gate_failed_groups   (문자열 리스트)
```

#### 진단 상태 확장

```
diagnosis_state ∈ {
    "scored",
    "insufficient_data",       # quality_score < 임계
    "conflicting_evidence",    # 모드 충돌
    "unsupported_chemistry",   # 화학 분기 미지원
    "insufficient_reference",  # 참조 코호트 부재
    "degenerate_fit",          # 파라미터 축퇴 (§5.6, §5.8)
}
```

**"모른다"를 1급 출력으로 만든다.** 점수를 강제로 내는 것보다 안전하다.

---

### 5.14 CC/CV 신호 기반 재검출

> **모듈:** `cyclediag/features/cc_cv.py` 재작성
> **버그:** raw에 `ChargeCVCapacity`가 있는데 feature는 `chgCVcapa=0` 보고.
> Si/ASSB에서 CV 구간은 **음극 수용 한계의 핵심 프록시**라 손실이 크다.

#### 알고리즘 — 라벨 비의존

```
1. 충전 leg 내에서 CV 구간 판정
   조건 (모두 만족):
     (a) V ≥ V_cutoff - 15 mV            # 전압 도달
     (b) dI/dt < 0                        # 전류 감쇠 중
     (c) |dV/dt| < 0.2 mV/min             # 전압 평탄
   → 연속 구간 중 가장 긴 것을 CV로 채택

2. 폴백
   raw에 ChargeCVCapacity 컬럼이 있으면 그 값과 대조하여 검증
   불일치 > 10 % → 경고 로그

3. 지표
   chgCVcapa    ← Q(CV 종료) - Q(CV 시작)
   chgCVtime    ← CV 지속 시간
   chgCCcapa    ← Q(CV 시작) - Q(충전 시작)
   chgCapa_CCratio ← chgCCcapa / (chgCCcapa + chgCVcapa)

4. CV 동역학
   I(t) in CV 구간에 I = I0·exp(-t/tau_CV) + I_inf 피팅
   tau_CV, I_inf_norm ← I_inf / I_cc

5. ★ 고정 CV 시간 기준 정규화 (중요)
   CV 종료 조건이 불명(#10/#11 UNKNOWN)하고 사이클마다 CV 시간이 변하면
   chgCapa_CCratio 는 열화가 아닌 프로토콜 변동을 반영한다.
   → 고정 기준 시간 T_ref (예: 30 min) 에서의 값으로 재정규화:
     Q_CV_at_Tref ← ∫_0^{T_ref} I(t) dt   (피팅된 지수함수로 외삽/절단)
     chgCapa_CCratio_norm ← chgCCcapa / (chgCCcapa + Q_CV_at_Tref)
   → 이 정규화 버전을 진단 입력으로 사용한다
```

#### 출력 컬럼

```
chgCVcapa, chgCVtime, chgCCcapa, chgCapa_CCratio,
chgCapa_CCratio_norm, Q_CV_at_Tref, tau_CV, I_inf_norm,
cv_detect_method ∈ {"signal", "column", "failed"},
cv_detect_mismatch_pct
```

#### 해석

동일 CC 조건에서 `tau_CV` 증가 + `chgCapa_CCratio_norm` 하락
= **음극 수용 능력 저하.** ASSB에서는 접촉 손실 또는 계면 저항 증가의 결과.
방전 방향 DC-IR만 있는 현 상황에서 **충전 방향 동역학의 최선 프록시**다.

#### 충전 방향 정보의 다른 우회로

```
chgR_asym_growth ← [ΔEoD_chgR / EoD_chgR_0] / [ΔR_DCIR_20 / R_DCIR_20_0]
```
`EoD_chgR_*`는 방전 후 rest 다음 충전 시작 구간의 ΔV/ΔI로, 충전 방향이고 매 사이클 나온다.

> ⚠️ SOC 불일치(EoD rest ≈ SOC 0 vs DC-IR SOC 20)이므로
> **절대 비율이 아니라 정규화된 성장률 비교로만** 사용한다.
> 이 비가 1을 넘어 커지면 충전 방향이 더 빠르게 열화한다는 신호.

**최소 프로토콜 추가 권장:** SOC 50에서 충전 펄스 30 s 1발 (1분 미만)
→ `ASR_asym = R_charge(50) / R_discharge(50)` 직접 지표 확보.

---

## 6. 파일 IO 견고화

### 6.1 step_type 라벨 의존 탈피

```
1차 (신호 기반, 주 판정):
   |I| ≤ rest_current_max            → rest
   I > 0 & |I| > rest_current_max    → charge
   I < 0 & |I| > rest_current_max    → discharge
   CV: 위 5.14의 3조건

2차 (라벨, 검증용):
   step_type과 불일치하면 경고 로그 + 사이클 flag
   → step_type을 required에서 optional로 강등
```

`rest_current_max` 기본값을 **셀 용량 비례**로 자동 설정: `0.007 × Q_rated` (72 Ah → 0.5 A).
현행 고정 0.01은 대형 셀에서 부적절하다.

### 6.2 프로토콜 자동 구조 인식

```
사이클별 특징 벡터:
   [mean_|I| / Q_rated, V_range, duration_h, has_rest, has_pulse, n_pulses]

클러스터링 (KMeans 또는 규칙 기반):
   → cycle_role ∈ {formation, routine_0.5C, RPT_C3, DCIR, unknown}

본 데이터에서의 기대 결과:
   |I| ≈ 25.78 A (0.33C) + 긴 duration        → RPT_C3
   |I| ≈ 38.67 A (0.5C)                       → routine
   |I| ≈ 77.34 A (1C) + 30 s + rest 3600 s    → DCIR
   TC1 rest only                              → formation/init

파생:
   protocol_block_id     연속 동일 role 블록
   baseline_cycle_auto   formation 이후 첫 RPT_C3
   rpt_anchor_cycles     RPT 블록의 cycle 2 목록
```

현행 `cycle_protocol.py`의 수동 제외 창을 대체한다.
`POST_RPT_EXCLUDE=5`는 유지하되 `protocol_block_id` 기반으로 자동 적용.

### 6.3 로더 견고화 체크리스트

- [ ] 인코딩 자동 감지 (cp949 / utf-8-sig / utf-16)
- [ ] 구분자, 멀티행 헤더, 단위 행 자동 감지
- [ ] **헤더에서 단위 파싱** — `Capacity(mAh)` → 스케일 자동
      ⚠️ 현행 "값이 크면 mAh" 휴리스틱은 72 Ah 셀에서 오작동 위험. **즉시 제거**
- [ ] 시간 단조성 / 중복 / 리셋 검사, resume 이음새 감지
- [ ] 결측 · 스파이크 필터
- [ ] `Temp (Celsius)` 전부 0.0 → `temperature_available=False` 로 명시 처리
- [ ] 채널 스왑 / 셀 ID 오기입 탐지 (용량 스케일 이상치)

### 6.4 성능

Ch22 raw ~80 MB. 폴더 배치가 병목이 될 수 있다.

```
캐시 키: sha256(파일) + code_version + config_hash
저장:    parquet (features), npz (dQ/dV 곡선)
병렬:    사이클 단위 (extract), 파일 단위 (batch report)
청크:    raw 로드 시 pandas chunksize + 사이클 경계로 재조립
```

---

## 7. 아키텍처 개선

### 7.1 화학 조성별 분기

현재 `config/mode_weights_fullcell_v1.json` **하나**로 모든 셀을 처리한다.
대상이 ASSB로 확인된 이상 물리적으로 성립하지 않는다.

```
config/
  mode_weights_assb_si_v1.json      ← 신규 (본 셀, §3.2 순위)
  mode_weights_liion_nmc_gr_v1.json ← 기존을 이쪽으로 이동
  chemistry_registry.json           ← 셀 메타 → config 매핑
```

셀 메타데이터 스키마:
```
cell_meta = {
  "electrolyte": "solid" | "liquid",
  "cathode": "S83S", "cathode_spec_cap_mAh_g": 215,
  "anode": "SJ-ASG903-1300", "anode_spec_cap_mAh_g": 1306,
  "anode_type": "Si_rich" | "graphite" | "LTO",
  "np_ratio": 1.08, "loading_cat": 5.0, "loading_ano": 5.4,
  "format": "?", "q_rated_Ah": 72,
  "stack_pressure_MPa": null,          # ← 최우선 확보 대상
  "chamber_temp_C": 45,
}
```

### 7.2 불확실성 전파

모든 지표가 **값 + σ**를 출력한다.

```
SoHQ            = 94.2 ± 0.3 %
R_ohmic         = 1.24 ± 0.09 mΩ      ← §5.3 피팅 공분산에서
LAM_curve_proxy = 6.1 ± 2.4 %         ← §5.6 자코비안에서
```

전파 경로:
```
v_noise_sigma (§5.13)
  → 피크 높이 σ (§5.8)
  → 모드 점수 σ (pattern_scoring)
  → confidence
```
σ가 큰 지표는 모드 점수에서 자동 가중치 하향.

### 7.3 참조 코호트

현행 "입력 배치 자체를 reference"는 순환 논리다.
동일 셀의 열화 진행분이 참조 분포에 섞여 anomaly를 희석한다.

```
cohort_key = (chemistry, format, protocol_hash, chamber_temp)
cohort_store/
  {key}/reference_stats.parquet    median, std per feature
  {key}/n_cells, n_cycles, version

reference 미존재 → diagnosis_state = "insufficient_reference"
                   (점수를 내지 않는다)
```

현재 Ch22/Ch25 2셀뿐이므로 **당분간 anomaly 점수는 참고용으로만** 표기한다.

### 7.4 반증 실험 추천 엔진

각 모드 점수에 "이 가설을 확정/기각하려면 무엇을 측정해야 하는가"를 첨부한다.

| 모드 | 추천 확인 |
|---|---|
| 접촉 손실 | **구속 압력 재인가 후 용량·R 회복 여부** (결정적), 단면 SEM |
| 계면 저항 증가 | EIS (Nyquist 반원 분리), 다중 펄스 길이 DC-IR |
| SE 분해 | XPS/XRD 계면 분석, 가스 발생량 |
| 관통/미세단락 | **장시간 rest 자가방전** (§5.4로 이미 산출), 절연 저항 측정 |
| LAM_PE | 반쪽셀 재조립, C/20 pseudo-OCV |
| LLI vs LAM | C/20 pseudo-OCV + 반쪽셀 템플릿 fitting |
| 고상 확산 제한 | GITT, rate capability sweep |

이 필드가 있으면 `cyclediag`가 "점수 출력기"에서 **"실험 설계 조력자"**로 격상된다.

### 7.5 인간 라벨 루프

golden set이 없는 현 상황의 유일한 탈출구.
현행 `pne_studio2` Diagnosis UI는 "displays only"다.

```
최소 구현:
  1. 진단 결과 옆에 [확인] [수정] [불확실] 버튼
  2. 수정 시: 실제 모드 + 근거 + 확인 방법 입력
  3. labels/{cell_id}_{cycle}.json 로 저장
  4. golden_set 빌드 스크립트가 수집

수십 케이스만 모여도 §7.1 가중치를 데이터 기반으로 보정 가능.
P0에서 착수해야 P5에서 쓸 수 있다. 가장 늦게 결실 맺는 투자.
```

### 7.6 재현성

```
모든 출력 행에 기록:
  code_version, config_hash, input_file_sha256,
  feature_set_version, random_seed

결과 diff 도구:
  두 실행의 feature 테이블을 비교해 어떤 컬럼이 얼마나 바뀌었는지 리포트
  → 코드 변경의 영향 범위를 자동 감지
```

### 7.7 이상 탐지의 시간 구조화

현행 z-score는 시간 무관 point anomaly만 잡는다.

```
추가:
  change_point   ← Bacon-Watts / PELT (급변 시점)
  trend_anomaly  ← rolling slope의 z-score
  collective     ← 여러 지표 동시 이탈 (Mahalanobis 거리)
```

### 7.8 BatteryML형 데이터·파이프라인 계층 (요약)

현행 `extract_features()` 평탄 테이블을 BatteryML식 계층으로 분리한다.
상세 로드맵·우선순위는 **[§12](#12-외부-오픈소스-참고-batteryml--pybamm)**.

```
CellData → CycleData[] → FeatureSet → Model → Evaluation(train/val/test)
```

원칙: **파서 / 공통 모델 / feature extractor / model / eval** 을 서로 독립 모듈로 유지.
YAML(또는 JSON) 실험 설정으로 재현성(§7.6)과 연결.

---

## 8. 검증 체계

golden set이 없으므로 **합성 데이터가 1차 검증 수단**이다.

### 8.1 합성 데이터 역복원 (즉시 구현)

```
forward model:
  1. 기준 곡선 V_ref(Q) 를 실측 BOL에서 채택
  2. 알려진 열화 적용:
       LAM: Q → s·Q          (s = 0.95)
       LLI: Q → Q + o         (o = 2 Ah)
       R:   V → V - I·dR      (dR = 5 mΩ)
       피크: 위치 shift + broadening + skew
       R(t): 알려진 (R_Ω, R_ct, tau_ct, A) 로 펄스 합성
       자가방전: 알려진 leak 전류로 rest 곡선 합성
  3. 노이즈 주입: v_noise_sigma 실측값, 양자화 적용
  4. 파이프라인 통과 → 복원 오차 측정

목표 복원 오차:
  s, o        < 2 %
  dR          < 5 %
  R_Ω, R_ct   < 5 %
  피크 위치    < noise_floor
  self_discharge < 10 %

→ 이것이 P2 이후 모든 변경의 회귀 테스트 기준선이 된다
```

### 8.2 공개 데이터셋 벤치마크

| 데이터셋 | 검증 항목 | 주의 |
|---|---|---|
| Severson (LFP 124셀) | `dQV_log_var` 재현, knee 예측, RUL | LFP이므로 dQ/dV 경로는 검증 불가 |
| Oxford Path | 곡선 정합, LLI/LAM 분해 (저율 OCV 포함) | 액체 셀 |
| Sandia (Preger) | 공변량 보정, 온도 의존성 | 액체 셀 |

**한계:** 공개 ASSB 장기 사이클 데이터셋은 희소하다.
→ 알고리즘의 **수치적 정확성**은 공개셋으로, **ASSB 해석 규칙**은 사내 검증으로 나눈다.

### 8.3 내부 교차검증 (독립 경로 일치도)

동일 물리량을 서로 다른 경로로 산출해 일치를 확인한다.

| 물리량 | 경로 A | 경로 B | 일치 조건 |
|---|---|---|---|
| 분극 증가 | §5.3 `R_30s_total` | §5.6 `fit_dR` | 같은 방향, 상관 > 0.8 |
| 확산 제한 | §5.3 `A_diff` | §5.10 `PER` | 같은 방향 |
| 시정수 | §5.3 `tau_ct` | §5.5 `tau_r1` | 같은 자릿수 |
| 손실 위치 | §5.6 `residual_argmax_SOC` | §5.9 `eta_argmax_SOC` | ±15 SOC 이내 |
| LAM | §5.6 `LAM_curve_proxy` | §5.8 `group_area` 감소 | 같은 방향 |

**불일치는 버그이거나 물리적 발견이다.** 자동 리포트에 포함한다.

### 8.4 Leave-one-RPT-out (§5.8 전용)

```
for each anchor pair (n_i, n_{i+1}):
    n_i 앵커 + 그 사이 0.5C 궤적으로 n_{i+1} 피크 예측
    실측 n_{i+1} RPT와 비교
    error ← RMS(예측 V_k - 실측 V_k)
    유의성 판정: error vs noise_floor (RPT cyc1-cyc2 차이)
```

RPT 주기 105 사이클이므로 이 검증이 특히 중요하다.

---

## 9. 실행 계획

### 9.0 Full-cell 우선 (하프셀 불필요) — **현재 메인 트랙**

§0의 F1–F8. 하프셀·PyDMA DMA는 §9.3 이후.

| 순서 | 기능 | 작업 | § |
|---|---|---|---|
| A | F1+F2 | ICA/DVA 생성 품질 고정 · 피크 검출 파라미터 스윕 | 5.1, 5.2 |
| B | F3 | peak matching (assign + evolution) 안정화 | peak_tracking, 5.8 |
| C | F4+F5 | 위치 이동 · 면적 감소 궤적 → pattern 입력 | peak_trajectory |
| D | F7 | R(t) 3성분 · polarization / PER | 5.3, 5.10 |
| E | F6 | ΔQ(V) · baseline curve correlation | 5.7, 5.6 |
| F | F8 | knee · change-point (SoHQ / R / peak_V) | 5.12, 7.7 |
| G | — | `mode_weights_assb_si_v1` + supporting_features 배선 | 7.1 |

### 9.1 즉시 (추가 실험 없이, 데이터 이미 존재)

| # | 항목 | § | 근거 |
|---|---|---|---|
| 1 | **dQ/dV 필터 스윕 진단** | 5.1 | §5.8 필요 여부를 결정. 다른 모든 피크 작업의 선행 |
| 2 | **R(t) 3성분 분해** | 5.3 | ASSB 접촉 손실 vs 계면 화학 분리. 10 Hz 확보됨 |
| 3 | **자가방전율** | 5.4 | 미탐지 실패 모드. rest 3600 s 이미 존재 |
| 4 | **`Q_relax` 전 블록 집계** | 5.10 | 노이즈 하한 0.065 % 확정됨 |
| 5 | **`chgCVcapa` 버그 수정** | 5.14 | Si 음극 수용 한계 프록시 복구 |
| 6 | **SOC 분해 히스테리시스** | 5.11 | Si 열화 직접 지표 |
| 7 | `R_ratio_20_50`, `R_SOC_slope` | 5.9 | 제한 전극 식별 |
| 8 | `VE`, `CI_per_hour` | 4.3 | 분극/쿨롱 분리, 45 °C 캘린더 성분 |
| 9 | 데이터 품질 지표 산출 | 5.13 | 이후 모든 진단의 전제. ADC 분해능도 역추정 |
| 10 | baseline → formation 후 첫 RPT | 6.2 | 현행 cycle=1은 구조적 편향 |

### 9.2 단기

| # | 항목 | § |
|---|---|---|
| 11 | 적응형 스무딩 / V축 보간 전환 | 5.2 |
| 12 | `R_recovery_tau` + `V_inf_est` | 5.5 |
| 13 | 3-파라미터 곡선 정합 | 5.6 |
| 14 | ΔQ(V) 통계 | 5.7 |
| 15 | η(SOC) × DC-IR 결합 | 5.9 |
| 16 | RCF, PER | 5.10 |
| 17 | 페이드 지수 · knee | 5.12 |
| 18 | 프로토콜 자동 구조 인식 | 6.2 |
| 19 | `mode_weights_assb_si_v1.json` 작성 | 7.1 |
| 20 | 합성 데이터 검증 하네스 | 8.1 |
| 21 | Q-domain 제약 디컨볼루션 | 5.8 **(§5.1이 (A)일 때만)** |

### 9.3 중기

| # | 항목 |
|---|---|
| 22 | 불확실성 전파 전면 도입 (§7.2) |
| 23 | 참조 코호트 관리 (§7.3) |
| 24 | 반증 실험 추천 엔진 (§7.4) |
| 25 | 인간 라벨 루프 (§7.5) |
| 26 | 재현성 인프라 (§7.6) |
| 27 | 공개셋 벤치마크 (§8.2) |
| 28 | 반쪽셀 확보 → forward/inverse fit → 실제 `*_est` |
| 29 | **BatteryML형 CellData→CycleData→FeatureSet 계층** (§7.8, §12.1) |
| 30 | **YAML 실험 설정 + train/val/test 평가 파이프라인** (§12.1) |
| 31 | **ΔQ(V) → RUL 예측 경로** (§5.7 + §12.1, Severson 재현 포함) |
| 32 | **전극 SOH 상태 벡터** (Q_PE / Q_NE / n_Li) — 스키마만, PyBaMM 비내장 (§12.2) |
| 33 | **OCP / CompositeOCP + fit_target(OCV/ICA/DVA)** — PyProBE API (§12.4) |
| 34 | **Batch DMA + warm-start + quantify_degradation_modes** — PyDMA 과학 (§12.3–12.4) |
| 35 | **stoich window · utilization · blend phase** 출력 스키마 (§12.3) |
| 36 | **ICA peak descriptor 스키마 정렬** — DiffCapAnalyzer (§12.5) |

> **철칙: 반쪽셀 또는 검증된 템플릿 확보 전까지 `*_est` / `*_est_hc_calibrated` 컬럼을 채우지 않는다.**
> null placeholder가 근거 없는 숫자보다 안전하다.
> **PyBaMM / PyDMA / PyProBE를 cyclediag 필수 의존성으로 내장하지 않는다** — 개념·API·알고리즘만 참고 (§12).

### 9.4 프로토콜 변경 제안 (외부 협의 필요)

| 항목 | 얻는 것 | 비용 |
|---|---|---|
| **충전 펄스 SOC 50, 30 s 1발** | `ASR_asym` → plating/관통 위험 직접 지표 | **1분 미만. 가성비 1순위** |
| **0.5C 펄스 1회 병행** | 1C 비선형성 검증 (§5.3) | 1분 |
| **온도 로그 export 복구** | Arrhenius, DTV, 셀 간 비교 | 설정 변경 |
| **구속 압력 기록** | `contact_loss_score` 해석의 전제 | 계측 추가 |
| 중간 경량 RPT (C/3 1사이클) | 앵커 간격 105 → 50 단축 | 사이클당 ~4 h |
| C/20 pseudo-OCV (RPT 5회당 1회) | 열역학 절대 기준, 반쪽셀 대체 부분 가능 | 고비용 |

---

## 10. 미해결 질문

### 10.1 최우선 (해석의 전제)

1. **구속 압력** — 초기값, 유지 방식(정압/정변위), 사이클 중 로깅 여부
   → ASSB 접촉 손실의 지배 변수. 없으면 `contact_loss_score` 해석 불가
2. **온도 실측** — export에서 `Temp` 복구 가능한지, 측정 위치(표면/분위기), 챔버 제어 밴드
3. **진단 결과로 내리는 실제 의사결정** — QC 스크리닝 / 기구 규명 / RUL 중 무엇인가
   → 미정이면 §9의 우선순위가 전부 동률이 된다. **기술이 아닌 제품 결정**

### 10.2 알고리즘 파라미터 확정에 필요

4. CV 종료 전류 (C-rate 또는 A) — §5.14 정규화 기준
5. `.sch` 샘플링 설정 원문 (Δt/ΔV/ΔQ 트리거) — §5.13 기대값 설정
6. ADC 전압 분해능 — §5.13에서 역추정 가능하나 확인이 나음
7. 셀 포맷 (파우치/각형) 및 정격 용량 — 설계 41.75/83.49 Ah vs 실측 72 Ah 불일치 원인
8. 전극 두께, Si 질량 %, 양극 조성비 — §7.1 메타 완성

### 10.3 검증 체계 설계에 필요

9. 코호트 규모 · replicate 수 — §7.3 참조 코호트 구성 가부
10. 의도적 열화 셀(저온 급속충전 등) 확보 가능성 — 라벨 케이스
11. 반쪽셀 제작 또는 문헌 OCV 템플릿 확보 가능성 — P3 도달 가부
12. 프로토콜 추가 시간 예산 — §9.4 실현 가능성
13. 현재 가장 신뢰받지 못하는 진단 출력 — 개선 우선순위의 직접 근거

---

## 11. 금지 사항

### DO NOT

- **액체 전해질 전제 모드**(전해액 고갈, wetting, 액상 확산)를 ASSB에 적용
- **`Q_NE_apparent`(흑연 stage)** 를 Si-rich 음극에 적용
- 피크 지표로 **LAM_NE 판정** (관측 피크는 양극 유래)
- **충전 곡선과 방전 곡선 직접 비교** (Si 히스테리시스, 경로 의존 OCV)
- 0.5C에서 뭉친 피크를 **강제로 8개로 분리**하고 개별 값 보고
- `R_DCIR`(30 s)을 **0.5C 정상상태 분극과 동일시** — PER > 1이 정상
- **1C 펄스 R 절대값**을 선형 응답 가정 하에 해석 (검증 전까지 상대 추이만)
- `R_10s/R_30s`를 **EIS 저항**이라 호칭
- `peak1`이 영구히 동일한 전기화학 반응이라 가정
- 데이터 품질이 나빠도 **점수를 강제 출력**
- 반쪽셀 확보 전 **`*_est` 컬럼 채우기**
- 단일 `mode_weights` config를 **모든 화학**에 적용
- **온도 로그가 비어 있는 상태에서 Arrhenius 보정** 적용
- RPT 앵커로 **cycle 1** 사용 (직전 0.5C 이력 보유 → cycle 2를 사용)
- `vp_diag`가 여전히 존재한다고 가정
- **BatteryML / PyBaMM / PyDMA / PyProBE / DiffCapAnalyzer를 필수 의존성·submodule로 내장** (개념·API만 참고, §12)

### DO

- 분리 불가를 **명시적 상태**로 보고 (`peak_k_resolved = False`, 그룹 단위)
- 보정 파라미터(α, s, o, ΔR, `Reff_scale`)를 **그 자체로 진단 지표**로 활용
- 모든 모드 점수에 **반증 실험** 첨부
- **RPT 2사이클 차이(0.065 %)를 노이즈 하한**으로 사용
- Layer 1(관측) → Layer 2(상태) → Layer 3(기구)를 **분리 표기**
- 독립 경로 **교차검증**(§8.3) 결과를 신뢰도에 반영
- 모든 결론에 **"45 °C 조건 한정"** 명시 (온도 로그 복구 전까지)
- `R_ohmic` 성장을 **접촉 손실의 직접 증거**로 해석 (ASSB 한정)


---

## 12. 외부 오픈소스 참고 (BatteryML · PyBaMM · PyDMA · PyProBE · DiffCapAnalyzer)

> **목적:** VP 진단 → 수명 예측 · **DMA(LLI/LAM)** · 전극 SOH로 확장할 때 **재발명하지 말고 구조를 빌린다.**
> **원칙:** 코드/패키지 통째 내장이 아니라 **과학 기능 · API 계약 · feature 정의 · 평가 파이프라인**만 흡수.
> **추천 조합:** **PyDMA의 과학적 기능** + **PyProBE의 API 구조** (§12.3–12.4). 수명 예측은 BatteryML (§12.1).
> **관련:** §5.6–5.8, §7.8, §8.2, §9.3 #29–36, [LLI_LAM_DIAGNOSIS.md](LLI_LAM_DIAGNOSIS.md) DM-P3

### 12.1 BatteryML에서 가져올 것

출처: [microsoft/BatteryML](https://github.com/microsoft/BatteryML)

| 아이디어 | cyclediag 적용 | 우선순위 | 연계 |
|---|---|---|---|
| **Cycler별 데이터 파서** | `io/` 에 vendor adapter (`pne`, 향후 `arbin`/`biologic`…) → 공통 스키마로 정규화. §6.3 로더 견고화와 동일 축 | 중기 | §6 |
| **공통 BatteryData 모델** | `CellData` → `CycleData[]` → (옵션) `StepData`. 현행 flat `extract_features()` 결과의 **상위 컨테이너** | **단기~중기** | §7.8, DATA_SCHEMA |
| **YAML 기반 실험 설정** | train/feature/model/split/metric을 한 파일로 고정. `config_hash`(§7.6)와 연결 | 중기 | §7.6 |
| **Feature extractor ↔ Model 분리** | `features/*` 는 벡터만, `models/*` 는 학습·추론만. diagnosis engine도 feature 입력만 받도록 | **단기** (이미 방향은 맞음, 계약 명문화) | ROADMAP Phase 1–2 |
| **train / validation / test 평가 파이프라인** | 셀 단위 hold-out (GroupKFold by `cell_id`). RUL·진단 공통 평가 CLI | 중기 | ROADMAP 0.5, 2.4 |
| **Severson ΔQ(V) feature** | 공통 V축 보간 → ΔQ(V)의 var/skew/kurtosis/min/`log_var`. **수명 예측 확장의 1순위 feature** | **단기** (알고리즘은 §5.7) | §5.7, §8.2 |

#### 목표 데이터 계층

```
CellData
  meta: cell_id, chemistry, protocol_hash, q_rated, ...
  cycles: list[CycleData]
    CycleData
      cycle_index, role (routine/RPT/DCIR/...), V/I/Q/t arrays (또는 lazy path)
      features: FeatureSet          ← extract_* 결과
FeatureSet
  version, family tags (§4.1), values + sigma (§7.2)
ModelSpec (YAML)
  features: [dQV_log_var, ...]
  model: ridge | xgb | ...
  split: {by: cell_id, train/val/test ratios}
EvaluationReport
  metrics, per-cell residuals, config_hash
```

#### 수명 예측으로의 확장 경로

```
1. §5.7 dQV_* 구현 + BatteryML Severson 공개셋 재현 (§8.2)
2. FeatureSet에 early-cycle window feature 묶음 (cycle 10–100 등)
3. YAML experiment → RUL / cycle-to-knee 회귀
4. ASSB(Ch22/Ch25)에 전이 — 화학·rate 차이로 계수 재학습, feature 정의는 공유
```

진단(LLI/LAM/접촉 손실)과 예측(RUL)은 **같은 FeatureSet**을 쓰되 **ModelSpec만 분기**.

### 12.2 PyBaMM에서 가져올 것

출처: [pybamm-team/PyBaMM](https://github.com/pybamm-team/PyBaMM)

| 아이디어 | cyclediag 적용 | 우선순위 | 비고 |
|---|---|---|---|
| **전극 SOH 상태 벡터** | 셀 전체 `SoHQ` 외에 `Q_PE_eff`, `Q_NE_eff`, `n_Li` (또는 LLI) 스키마. Layer 2 출력과 정렬 | 중기 (스키마) / 장기 (추정) | electrode SOH solver 예제 참고 |
| **열화 모드 온톨로지** | SEI · plating · LAM · LLI를 **파라미터 이름 공간**으로 문서화. ASSB에서는 §3에 맞게 재매핑 (액상 SEI≠고체 계면상) | 문서·스키마 | §3, LLI_LAM_DIAGNOSIS |
| **Forward model 합성** | §8.1 합성 데이터 검증의 물리 근거. 필요 시 **오프라인**으로 PyBaMM 곡선 생성 → cyclediag regression fixture | 선택 | **런타임 의존성 금지** |
| **Half-cell / OCV 라이브러리 사고방식** | DM-P3 half-cell calibrate와 동일 철학: 템플릿 OCV → full-cell 정합 | 장기 | §5.6, DM-P3 |

#### 명시적 비목표

```
✗ PyBaMM을 cyclediag requirements에 추가
✗ SPM/DFN을 진단 루프 안에 돌림 (너무 무겁고 ASSB 파라미터 부재)
✗ 액체 전해질 기본 파라미터를 ASSB에 그대로 적용
```

가져올 것은 **"셀 용량 하나"가 아니라 전극·리튬 재고 상태를 분리해 보고하는 출력 계약**이다.
수치 해는 반쪽셀·C/20·곡선 정합(§5.6) 경로로 채우고, PyBaMM은 **검증용 forward 시뮬레이터(옵션 툴)** 로만 둔다.

### 12.3 PyDMA에서 가져올 것 (DMA 과학 — **목표에 가장 근접**)

출처: [tum-ees/PyDMA](https://github.com/tum-ees/PyDMA)

full-cell **pseudo-OCV**와 양극·음극 **half-cell OCP**를 맞춰 다음을 **정량화**하는 구조가 cyclediag DM-P3 / Level 3과 동일한 문제 정의다.

| 정량 출력 | cyclediag 매핑 | 비고 |
|---|---|---|
| **LLI** | `LLI_est_hc_calibrated` (Level 3) | 하프셀·저율 OCV 확보 후 |
| **LAM at anode / cathode** | `LAM_NE_est_*`, `LAM_PE_est_*` | ASSB Si-rich: LAM_NE는 피크 단독 금지 (§3) — **OCP fit 경로로만** |
| **전극 utilization · stoichiometry window** | `stoich_window_PE/NE`, `utilization_*` | 최근 PyDMA: blend phase별 window 확인 |
| **전극 불균일성** | `inhomogeneity_*` (신규 지표 후보) | pattern score 보조 증거 |
| **Si–graphite blend 변화** | `CompositeOCP` / blend fraction drift | SJ-ASG903 Si-rich blend에 **직접 관련** |

#### 과학 기능 (흡수 대상)

```
1. OCV뿐 아니라 DVA(dV/dQ) · ICA(dQ/dV)를 가중치로 동시 fitting
2. 양극·음극 모두 blend electrode로 모델링 가능
3. 각 전극 / blend phase의 stoichiometry window 검증
4. pseudo-OCV(full-cell) ↔ half-cell OCP 라이브러리 정합
```

#### ASSB 적용 시 주의

- 현재 SJ900은 **반쪽셀 OCV 없음** → Level 3 수치는 placeholder. PyDMA식 fit은 **템플릿 OCP 확보 후** 활성화.
- C/3 RPT를 pseudo-OCV 근사로 쓸 때는 §5.9·§5.13 품질 게이팅 필수 (분극 잔류).
- Si 히스테리시스: **충·방전 각각** fit하거나, PyProBE식 평균 OCV로 저항 보정 후 fit (§12.4).
- 액체 셀 기본 OCP를 ASSB SE 계면에 그대로 쓰지 말 것.

**우선순위:** 하프셀 또는 문헌 OCP 템플릿 확보 시 **DM-P3 1순위 참고 구현**. 그 전에는 API·출력 스키마만 정렬.

### 12.4 PyProBE에서 가져올 것 (API · 배치 DMA 구조)

출처: [ImperialCollegeLondon/PyProBE](https://github.com/ImperialCollegeLondon/PyProBE)

실제 DMA 모듈이 있고, half-cell OCP → full-cell fitting으로 SOH / LAM_PE / LAM_NE / LLI를 낸다.
**과학은 PyDMA, 모듈 경계·배치 UX는 PyProBE** 를 기본 추천으로 둔다.

| API / 기능 | cyclediag 적용 |
|---|---|
| **`OCP` / `CompositeOCP` 객체** | `diagnosis/halfcell/` — 단일·blend OCP 타입. Si–Gr blend = CompositeOCP |
| **OCP data interpolation** | 공통 V 또는 Q 그리드로 half/full 정렬 (§5.6·§5.7과 동일 패턴) |
| **Fitting target 선택** | `fit_target ∈ {OCV, ICA(dQ/dV), DVA(dV/dQ), weighted}` — config로 전환 |
| **Batch DMA (여러 RPT)** | `rpt_anchor` 주기마다 fit; 병렬·순차 스위치는 기존 적응형 병렬 정책과 공유 |
| **이전 RPT → 다음 fit 초기값** | warm-start: `x0_{k} ← x*_{k-1}` — 105사이클 앵커 간격에서 수렴 안정화 |
| **충·방전 OCV 평균** | 저항(분극) 영향 보정 후 열역학 곡선에 가깝게 — Si 히스테리시스와 **병기** (평균만으로 히스테리시스를 지우지 않음) |
| **`quantify_degradation_modes()`** | 여러 OCV fit 결과 → `LLI`, `LAM_PE`, `LAM_NE` 직접 반환. Level 2/3 export 계약의 참고 시그니처 |

#### 목표 모듈 스케치

```
diagnosis/halfcell/
  ocp.py              # OCP, CompositeOCP (PyProBE 스타일)
  interpolate.py
  fit_ocv.py          # target: ocv | ica | dva | weighted (PyDMA 가중 아이디어)
  dma_batch.py        # multi-RPT, warm-start, parallel/sequential
  quantify.py         # quantify_degradation_modes(fits) → LLI, LAM_PE, LAM_NE, windows
```

Level 1 pattern score와 **병행**: half-cell 없이도 pattern은 유지하고, DMA fit은 `diagnosis_state=halfcell_ready`일 때만 `*_est_hc_calibrated`를 채운다.

### 12.5 DiffCapAnalyzer에서 가져올 것 (ICA 피크 descriptor)

출처: [nicolet5/DiffCapAnalyzer](https://github.com/nicolet5/DiffCapAnalyzer)

오래된 프로젝트이나 **differential capacity 피크 descriptor** 구조는 peak tracking / §5.8과 맞는다.

| descriptor | cyclediag 컬럼 후보 | 연계 |
|---|---|---|
| peak voltage | `peak_*_V` | 기존 |
| peak height | `peak_*_H` / `H_norm` | PEAK_TRACKING |
| peak area | `peak_*_area` / `group_area` | §5.8 |
| peak width | `peak_*_W` | 기존 |
| positive / negative peak | `peak_*_sign` (ICA 부호) | 충·방전 분리 저장 |
| cycle별 피크 변화 | trajectory / delta 테이블 | peak_trajectory |

추가 참고:
- **Gaussian baseline fitting** — 검수 UI에서 baseline on/off (pne_studio Diagnosis / peak review)
- **cycle별 결과 저장** — golden·review parquet 스키마에 descriptor 고정 길이 벡터

**우선순위:** 단기 — 기존 `dqdv_peaks` / review export의 **스키마·네이밍 정렬**. 새 피크 엔진 교체는 아님.

### 12.6 구현 체크리스트 (로드맵 항목)

**데이터·수명 (BatteryML / PyBaMM)**
- [ ] `schema/cell_cycle.py` — `CellData` / `CycleData` / `FeatureSet` dataclass
- [ ] `io/adapters/` — PNE 파서 → CycleData
- [ ] `features/extract.py` — flat + FeatureSet 동시 출력
- [ ] `features/dqv_stats.py` — §5.7; Severson 교차검증
- [ ] `config/experiments/*.yaml` + `pipeline/evaluate.py`
- [ ] Layer 2 placeholder: `Q_PE_eff` / `Q_NE_eff` / `n_Li`
- [ ] `tools/synth_from_pybamm.py` (optional extra only)

**DMA (PyDMA 과학 + PyProBE API)**
- [ ] `diagnosis/halfcell/ocp.py` — `OCP`, `CompositeOCP` (+ blend phase window)
- [ ] OCP interpolation + `fit_target` config (OCV / ICA / DVA / weighted)
- [ ] `dma_batch.py` — multi-RPT, warm-start, parallel|sequential
- [ ] `quantify_degradation_modes()` — LLI, LAM_PE, LAM_NE, stoich windows, utilization
- [ ] 충·방전 평균 OCV 옵션 (저항 보정) vs 방향별 fit (Si 히스테리시스) 분기
- [ ] ASSB: Si–Gr CompositeOCP 템플릿 슬롯 (데이터 확보 전 null)

**ICA 피크 (DiffCapAnalyzer)**
- [ ] peak descriptor 스키마 고정: V, H, area, W, sign, cycle-delta
- [ ] review UI / export에 Gaussian baseline 토글 (선택)
- [ ] PEAK_TRACKING · §5.8 group metrics와 컬럼 정렬

### 12.7 하지 말 것 (이 절 한정)

- BatteryML / PyBaMM / PyDMA / PyProBE / DiffCapAnalyzer를 **git submodule·필수 의존성으로 벤더링**하지 말 것
- 공개 액체 셀 OCP·모델 가중치를 **ASSB 진단 점수에 직접 복사**하지 말 것
- ΔQ(V)만으로 ASSB 기구(접촉 손실 등)를 **단정**하지 말 것 — RUL 보조 feature; 기구는 §3·§5.3
- 반쪽셀·검증 템플릿 없이 PyDMA/PyProBE식 **`*_est_hc_calibrated`를 숫자로 채우지** 말 것
- 충·방전 OCV **평균만**으로 Si 히스테리시스를 “제거된 열역학”이라 단정하지 말 것
- DiffCapAnalyzer UI/스택을 **통째 이식**하지 말 것 — descriptor·baseline 개념만
