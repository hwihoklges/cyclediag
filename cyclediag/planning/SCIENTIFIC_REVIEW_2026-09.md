# 과학적 정확성 계약과 검증 로드맵 — 2026-09-18

대상: `review/correctness-and-roadmap`. [기존 로드맵](ROADMAP.md), [상세 개선안](IMPROVEMENT_ROADMAP.md), [peak 추적 계획](PEAK_TRACKING_ROADMAP.md)은 연구 이력으로 보존한다. 과거의 구현 완료 표시는 외부 과학 검증 완료와 다르다. 본 문서는 이번 리뷰의 해석 경계와 앞으로의 승인 기준을 명시한다.

## 1. 구현 상태와 해석 계약

**구현됨:** 단위 정규화·경고, 결측/상수 통계의 방어, derivative·IR 관련 수치 조건 검사, 품질 evidence coverage와 진단 validity metadata, proxy 점수의 휴리스틱 표기, 회귀 테스트. **실험적:** full-cell 메커니즘 추론, ECM 성분 해석, peak 화학 라벨 및 미교정 점수. **계획:** 독립 측정 기반 외부 검증, 교정된 불확실성, 결과 전 경로의 재현성 manifest. 테스트 성공은 알고리즘 계약의 증거이지 열화 원인의 증명이 아니다.

### 단위와 데이터 출처

- 내부 표준은 전류 **A**, 용량 **Ah**, 전압 **V**, 시간 **s**다. `f_Q_spec`은 active mass(g) 기준 **mAh/g**이며 Ah/g와 1,000배 차이가 난다. 질량의 대상·측정 방법도 기록해야 한다.
- `ColumnMap.units`는 논리 필드명의 입력 단위를 명시한다. 단위 없는 열을 명시적으로 설명할 수 있지만, 헤더의 단위와 충돌하면 `ValueError`로 거부한다. 차원이 다른 단위, 지원하지 않는 단위 표기 및 알 수 없는 논리 필드도 거부한다. capacity override는 개별 capacity 단위가 없을 때 Ah 계열 필드의 기본값으로 사용한다.
- 정규화 metadata가 이미 있는 열에 비표준 배율 override를 재적용하지 않는다. 단위 없는 입력은 수치 크기로 추측하지 않고 `assumed_canonical` 및 `unit_assumed` 경고를 기록한다. 이것은 단위가 검증되었다는 뜻이 아니다.
- `DataFrame.attrs`의 `unit_metadata`, `unit_warnings`, `unit_schema`, `unit_missing_fields`는 **그 자체로 영속 저장 계약이 아니다**. 일반 CSV, 재구성, join 등에서 소실될 수 있다. 이제 [checked feature CSV sidecar](REPRODUCIBILITY.md)는 명시적 serialization과 재로딩 무결성 검증을 제공한다. 이 helper 밖에서 attrs만으로 추적 가능하다고 주장하지 않는다.

### 지표별 해석 제한

| 지표 | 사용 계약 | 금지할 해석 / 추가 검증 |
|---|---|---|
| CE | 같은 cycle·동일 종료 조건의 방전량 / 전체 충전량 ×100(%). CC와 CV 충전량 포함 여부를 명시하고 이미 합계인 카운터에 CV를 이중 가산하지 않는다. | 100% 초과를 clipping하지 않는다. 측정 오차·카운터 reset·구간 불일치·프로토콜 문제를 조사하며 특정 메커니즘 증거로 단정하지 않는다. |
| VE | 현재 출력은 fraction(비율), percent가 아니다. 전압 평균 방식과 비교 구간을 기록한다. | VE 자체가 에너지 효율은 아니다. CE를 fraction으로 환산하고 같은 구간의 capacity-weighted 평균 전압을 쓴 경우에만 CE×VE가 에너지 효율과 대응한다. 별도로 V·I 시간 적분을 대조한다. |
| IR / DCIR | 기본은 ΔV/ΔI(Ω). 전류 부호, pre-pulse 기준, 측정 지연·분해능, SOC·온도·완화 상태를 일치시킨다. | 서로 다른 SOC/온도/휴지시간의 값을 직접 비교하지 않는다. ECM의 R_ohmic/R_ct/확산 항은 fit 성분이며 유일한 물리 메커니즘 분리가 아니다. |
| dQ/dV·dV/dQ | 단위, leg, 재표본화 간격, smoothing window·order, trim, peak prominence·matching config를 함께 보존한다. | 미분은 노이즈를 증폭한다. smoothing·보간·분해능·plateau가 peak를 생성/이동/소실시킬 수 있다. peak 하나를 LLI/LAM이나 상전이에 바로 대응시키지 않는다. |
| hysteresis·entropy·phase | 동일 SOC 정렬·온도·rate·relaxation 조건과 기준 전극 정보를 확인한다. | hysteresis에는 kinetics가 섞이며 단일 온도 곡선으로 entropy를 추정하지 않는다. 평형 전위의 온도 의존성·충분한 안정화가 없는 phase/entropy 해석은 실험적이다. |
| LLI/LAM/kinetics | full-cell의 중첩된 proxy를 설명 가능한 패턴 증거로 사용한다. | full-cell만으로 원인을 유일하게 식별하거나 절대 손실 %를 확정하지 않는다. 서로 상관된 feature는 독립 증거 여러 개가 아니다. |

`unknown`은 건강함이 아니다. 측정 부재, reference 부족, zero variance, domain 불일치와 실제 양호를 구분한다. `pattern_score`·confidence는 **교정되지 않은 heuristic support**이며 고장 확률·메커니즘 확률·신뢰구간으로 읽으면 안 된다. 현재 validity는 `causally_identified=false`, `validated_probability=false`로 제한을 알린다.

## 2. 현재 품질 게이트의 실제 범위

[quality 구현](../features/quality.py)은 사용 가능한 항목의 clipped ratio를 가중 기하평균한다. 목표/가중치는 samples_per_mV 0.5/0.25, voltage span/noise proxy 10/0.25, voltage coverage 0.9/0.20, rest/τ 3/0.15, pulse 안정도 1/0.15이다. 결측 항목을 좋은 값으로 채우지 않고 사용한 가중치 비율을 `quality_evidence_coverage`로 기록한다.

- `quality_status`: 근거 없음 `unknown`, 관측 목표 미달 `failed`, 미달은 없지만 일부 근거만 있음 `partial`, 전체 항목 평가 `assessed`. 높은 부분 점수는 충분한 데이터나 진단 적합성을 보증하지 않는다.
- `leg_completeness`라는 이름에도 현재 계산은 cycle 전체 전압 span / 기대 창(default 2.5–4.2 V)이다. **충전·방전 각각의 완결성 검사가 아니다.** `dqdv_snr`도 실제 derivative SNR이 아니라 voltage span/noise proxy다.
- rest 충분성은 측정/제공된 양의 τ가 있어야 계산한다. 기본 rest 전류 0.5 A, pulse 검출 `max(rest_current_max×10, 20 A)` 등은 특정 데이터 스케일의 휴리스틱이다. 작은 셀이나 다른 protocol에 일반화하지 않는다. step 경계·counter reset을 엄격하게 분리하는 pulse 품질 판정은 추가 검증 대상이다.
- [진단 constraints](../diagnosis/constraints.py)는 실패 그룹이 있거나 quality score≤0이면 mode confidence를 0으로 줄인다. 제외 protocol도 confidence를 0으로 만든다. 온도 부재·압력 불명·half-cell 미교정은 별도 감점/경고다. **품질 결측 자체가 모든 진단 경로를 강제 차단하는 것은 아니다.** indicator Track A와 메커니즘 Track B의 처리가 같다고 가정하지 않는다.
- 따라서 다음 단계는 프로토콜별 필수 evidence를 선언하고 unsupported/unknown을 보고서까지 보존하는 것이다. 보편적 임계값 또는 검증된 합격 판정으로 현재 score를 홍보하지 않는다.

## 3. 독립 검증 설계와 단계별 인수 기준

| 단계 / 상태 | 선행 조건 | 산출물 | 인수 기준 |
|---|---|---|---|
| C0 데이터 계약 — 부분 구현, 영속화 계획 | 현재 unit normalization | versioned schema, 단위·cell/batch·cycle/step·protocol·온도·질량 출처, attrs sidecar roundtrip | A/mA, Ah/mAh, V/mV, s/min 변환 불변성; conflict 거부; 재정규화 무중복; 결측·counter reset·중복/역전 시간 fixture; 재로드 후 단위 근거 동일. |
| C1 검증된 feature 추출 — 부분 구현, 검증 계획 | C0 | feature별 구간 정의·config snapshot, 품질 필수 항목, peak/IR/CE/VE reference benchmark | 독립 적분 용량·전압 계산과 사전 정의 오차 한계 내 일치; 노이즈/smoothing 민감도 공개; 불충분 조건에서 unknown; cell·batch hold-out에서 성능 유지. |
| C2 외부 검증 진단 — 실험적 → 계획 | C1, 승인된 reference dataset | full-cell proxy와 독립 측정의 비교, domain card, calibrated uncertainty 별도 출력 | 용량·OCV·적절한 half-cell/EIS 증거와 연결; 원인 구별 실패 사례 공개; held-out cell 및 batch에서 coverage/calibration·오탐/누락 보고. 미교정 support와 확률/구간을 명확히 분리. |
| C3 재현 가능한 batch/cloud·보고서 — 일부 도구 구현, 통합 계획 | C0/C1, 진단 주장은 C2 필요 | versioned job/result manifest, 제한 저장소·quota, 재시도/취소, 결과 보고서 및 scheduler 교환 계약 | 동일 input/config/code로 허용 오차 내 동일 결과, 병렬/순차 동일성, 중단 재시작·cross-project 접근 방지, 모든 unknown/경고·제외 근거 보존. 실행 중 장비 파라미터 자동 변경 금지. |

### 승인 reference dataset 절차 (계획)

1. 데이터 관리자가 원본 사용 권한, instrument calibration, chemistry, batch/cell 식별자, protocol·온도·SOC·압력·측정 오차와 hash를 확인하고 versioned manifest를 승인한다. 기존 DOE/halfcell fixture의 존재만으로 승인 완료가 아니다.
2. 용량 reference는 독립 적분/교정 측정으로, OCV는 충분히 완화된 조건에서 얻는다. Half-cell은 BOL/aged 대응·전극 loading·조건 일치와 채취 영향까지 기록한다. BOL OCP만으로 aged 절대 손실을 입증하지 않는다. EIS는 SOC·온도·선형성·정상성·주파수 범위가 적절할 때만 비교하며 ECM label의 ground truth로 취급하지 않는다.
3. 셀 단위 및 batch 단위 hold-out을 고정한다. 같은 셀의 다른 cycle이 train/test 양쪽에 들어가면 leakage다. golden 선정, scaling, smoothing 튜닝, peak matching 학습, threshold와 calibration은 train/validation 내부에서만 한다. test는 최종 평가까지 봉인한다.
4. ground-truth synthetic는 변환·counter reset·결측·peak overlap·known parameter 회수 검증에 사용한다. **합성 데이터 성공은 실제 셀의 메커니즘 증명이 아니다.** 실제 측정 benchmark를 별도 유지한다.
5. 분석 전에 feature 허용 오차, 실패/abstention 조건, uncertainty coverage 목표와 표본 수 근거를 도메인 책임자가 승인한다. 검증되지 않은 보편 수치 목표나 문헌 DOI를 만들어 넣지 않는다. 인용할 문헌의 DOI·원문·실험 조건을 확인하는 bibliography review는 미완료 항목이다.

## 4. 재현성 상태와 CI

**현재 구현:** config dataclass/JSON, unit·validity metadata, CI dependency snapshot과 [bounded provenance/checked CSV sidecar](REPRODUCIBILITY.md). `extract_features`/`diagnose_csv` 및 지원 CLI 경로는 input SHA-256, 기본값·단위·selector 포함 config hash, code content hash/commit/dirty, 제한 dependency 버전을 기록한다. **남은 작업:** 모든 저수준·batch·peak·model 경로의 적용 config/환경/random seed/데이터·모델 버전/exclusion 목록을 통합한다. 일부 hash나 attrs가 존재한다고 end-to-end 재현성이 완성된 것은 아니다.

이번 bounded provenance 추가 후 로컬 전체 suite는 RuntimeWarning 오류 모드에서 **221 passed**, 새 provenance 회귀 테스트는 **12 passed**다. 이는 원격 matrix 또는 실제 LFS 검증 결과가 아니다.

최종 로컬 Python 3.14.3 환경에서 `reports` 의존성을 포함한 전체 테스트를 RuntimeWarning 오류 모드로 재실행하여 **221 passed / 0 skipped**를 확인했다. 원격 matrix 및 실제 LFS 통합 통과를 이 결과로 대신하지 않는다.

- [cyclediag-ci.yml](../../.github/workflows/cyclediag-ci.yml): Ubuntu/Python 3.10, Ubuntu/3.14, Windows/3.14에서 **전체** test directory를 수집한다. 6개 파일만 선택하던 구성을 제거했다. LFS checkout/pull 없음, `-rs` 및 JUnit으로 fixture skip을 공개한다. RuntimeWarning을 오류로 취급한다.
- `reports` extra는 python-pptx를 설치하여 portable slide 검증이 의존성 부재로 빠지지 않게 한다. Windows Origin/PowerPoint COM 자동화 검증과는 다르다. `origin`의 pywin32는 Windows에만 설치된다.
- [cyclediag-integration.yml](../../.github/workflows/cyclediag-integration.yml): 수동 `workflow_dispatch`에서 `run_integration=true`일 때만 M01Ch022 raw **한 파일**을 selective LFS pull한다. pointer 여부를 검사한 뒤 presentation PNG/CSV 생성을 확인한다. PR마다 전체 약 821 MB fixture를 받지 않는다.
- 모든 workflow는 읽기 전용 권한, 중복 취소, timeout, 14일 artifacts를 사용한다. pip freeze와 JUnit(단위 job), presentation 산출물을 보존한다. 원격 실행 결과는 아직 확인 대상이다.
- 검증된 환경의 [Python 3.14 constraints](../../constraints/py314.txt)를 **3.14 job에만** 연결했다. 버전 constraints와 freeze는 배포물 hash를 검증하는 lock의 대체물이 아니다. 최소 Python 버전과 최신 버전 모두의 원격 의존성 해석·테스트를 확인해야 한다.