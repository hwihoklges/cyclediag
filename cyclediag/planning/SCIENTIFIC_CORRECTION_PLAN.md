# 과학적 보정·검증 실행 계획 (2026-09-19)

기준 코드: `834a15f`. 상태: **독립 검토 조건 수용·구현 승인 (2026-09-19)**.
계산 오류 수정, 물리 모델 교정, 실제 원인 검증을 구분한다. 코드 테스트는
차원·계산·거부 조건을 검증하며 실제 셀의 열화 원인을 입증하지 않는다.

## 1. 방법과 확인된 문제

| ID | 확인된 문제 | 수정 방법 | 합격 기준 |
|---|---|---|---|
| C1 | 증가하는 방전 누적량을 충전량과 같은 SOC 방향으로 비교 | 충전 s=(Q-Q0)/ΔQ, 방전 s=1-(Q-Q0)/ΔQ를 명시적으로 정렬; 취득 순서에서 reset/NaN 거부 | 동일 가역 선형 V(SOC)에서 차이 0; 0.1 V 상수 차이에서 2–98% 영역 적분 0.096 V; 용량 크기·상수 offset 불변 |
| C2 | Ea=20 kJ/mol로 온도 보정을 자동 적용 | 기본 보정 비활성화, raw 값 보존; 명시적 Ea·기준 온도·온도 범위·calibration id·측정 온도 출처가 있어야 opt-in | 기준 온도에서 항등, 명시적 Ea=0 항등, 역변환 검증, 비유한/절대영도 이하/범위 밖 거부 |
| C3 | counter 최대값 용량과 관측 구간 에너지 혼합 | 기존 reported total은 호환용 보존; observed ΔQ·E와 계산 구간을 별도로 기록; 동일 구간 VE/CE/EE 관계 유지 | E=VΔQ, observed CE×VE=EE (CE는 fraction 환산), counter offset 불변; 불완전 cycle을 완전 cycle로 주장하지 않음 |
| C4 | 측정 품질·chemistry 미검증이어도 diagnosis_valid 가능 | 수치 support와 해석 eligibility 분리; chemistry 명시 일치, 필수 측정 그룹/프로토콜 충족 여부를 구조화; legacy 탐색 사용은 unverified | missing/partial/failed/NaN/문자열 boolean을 안전하게 처리; 불충분 근거에서 판정 보류·확률 주장 없음 |
| C5 | 고정 0.5 A/20 A와 step_time reset을 넘는 pulse 집계 | 명시적 절대 임계값 또는 정격용량×명시적 C-rate를 사용; 미설정은 legacy 가정 경고; 첫 연속 pulse와 rest 구간별 시간·step 경계 처리 | 서로 다른 크기의 셀에서 동일 C-rate 입력 불변; 복수 pulse를 합쳐 안정도 계산하지 않음; reset/결측 구간 분리 |

### 중요 해석 경계

- C1의 s는 **endpoint-normalized SOC proxy**다. 서로 다른 cutoff, 미완료 leg,
  불균형 용량, drift가 있으면 실제 같은 SOC라는 증거가 아니다.
  부하가 걸린 V 차이를 평형 thermodynamic hysteresis로 부르지 않는다.
- C2의 Arrhenius 변환은 측정 조건/SOC/time scale별 경험 모델이다.
  calibration id가 있다는 사실만으로 교정 실험을 검증했다고 하지 않는다.
  측정 온도의 위치·시점이 불명확하면 보정하지 않는다. Ea를 임의 추정하지 않는다.
- C3에서 nonzero 시작 counter가 단순 offset인지 초기 샘플 누락인지 코드만으로
  구별할 수 없다. observed-window 지표와 reported total의 의미를 분리한다.
- C4에서 high support도 인과 확률이나 LLI/LAM 절대량이 아니다. chemistry 충돌은
  계속 거부한다. 모든 프로토콜에 rest/pulse를 의무화하지 말고 필수 evidence를
  명시적인 configuration으로 선언한다. 없는 evidence를 통과로 채우지 않는다.

### C4 코드 상태 (2026-09-23)

- C4 **구현 및 synthetic 회귀 테스트 완료** (Windows: 281 passed, `-W error::RuntimeWarning`).
  mode별 JSON `eligibility.required_groups`는 필요한 evidence 대안과 관련 quality 지표를
  선언한다. 해당 지표의 결측/실패는 거부하지만 무관한 rest/pulse 실패는 전역 거부가 아니다.
- `diagnosis_valid`는 이제 heuristic support와 별도로 scientific eligibility를 통과한
  행/모드만 참이다. 기존 점수·confidence는 탐색용으로 보존하고
  `diagnosis_heuristic_support_valid`, 모드별 `*_scientific_eligible` 및
  `*_eligibility_reasons`로 구분한다. JSON sidecar 및 진단 CSV에도 결과를 기록한다.
  `scientific_validity.per_mode`에 상세 그룹/이유를 보존하며 contract version은
  `mode_eligibility_v1`이다. 기존 `fullcell_v1`/`assb_si_v1` 모델 가중치와 점수
  의미는 바꾸지 않았다; 구 CSV의 `diagnosis_valid`와 새 필드는 동일 의미가 아니다.
- protocol comparability (`protocol_comparable=True`, `protocol_kind=routine`,
  `protocol_excluded=False`, baseline 사용 시 동일한 명시적 `protocol_id`) 및
  chemistry 일치, 유효한 quality/feature/baseline/fit metadata가 있어야 eligibility가
  성립한다. 기존 추출 데이터에 comparability 증명이 없으면 탐색 점수는 유지하되
  **unverified**로 보류한다. 문헌/실측 calibration이나 실제 원인 확인은 아직 없으며
  eligibility는 검증된 인과 판정·확률·절대 열화량이 아니다.
- C5 정격용량 scaling 역시 보편적인 검증 임계값은 아니다. 사용자가 제시한
  current/C-rate 기준과 근거를 기록하며 센서 분해능·protocol별 검증이 필요하다.

### C5 코드 상태 (2026-10-01)

- 품질 지표의 pulse 기준은 명시적인 A 또는 명시적 정격용량(Ah) × C-rate(h⁻¹)로
  결정한다. 직접 호출에 기준이 없으면 pulse 지표는 null이며 이전 고정 20 A 기준은
  제거했다. 추출 경로의 기본 72 Ah·1C는 호환성을 위해 사용하지만
  `legacy_protocol_default_unverified` 출처·경고를 남긴다. 프로토콜별 임계값은
  아직 교정/승인되지 않았다.
- 첫 연속 pulse만 1초 샘플 수 및 전류 안정성에 사용한다. step 변경, 부호 변경,
  누락, 같은/역행 시간에서 구간을 분리하며 rest 역시 같은 경계를 적용한다.
  `step_time`이 없으면 초 단위 결과를 추정하지 않는다. 이는 측정된 pulse의
  전압 step 응답이나 true pre-pulse ΔV/ΔI 검증을 대신하지 않는다.
- synthetic 회귀에 임계값 차원, 작은 셀, 다중 pulse, step/clock reset, 시간 결측을
  포함한다. 문헌·센서·장비 데이터 교정 및 독립 성능 검증은 미완료다.

## 2. 검토에서 배제할 잘못된 보정

1. C-rate 시간은 t=3600f/C (s)이다. 동일한 nominal 기준에서 용량이 약분된다.
   50 mAh 셀도 0.5C로 nominal 전량을 처리하는 이상 CC 시간은 7200 s다.
2. CE>100%는 과충전 또는 특정 고장의 증거가 아니다. raw 측정값을 보존하고
   계측·구간·초기 상태를 조사한다. clipping하거나 causal evidence로 사용하지 않는다.
3. counter reset을 정렬·smoothing으로 숨겨 에너지를 만드는 보정은 금지한다.
4. chemistry mismatch를 임의 완화하거나 어떤 degradation mode가 보편적이라고
   지정하지 않는다. 기존 evidence coverage는 이미 있으며 별도 확률로 바꾸지 않는다.
5. 검증되지 않은 penalty·Ea·불확실성 범위를 새로 만들어 넣지 않는다.

## 3. 실행 순서와 검토 관문

1. 소스·실제 호출 경로·기존 테스트 대조 (**완료**).
2. 본 계획을 별도 검토자로 점검: 차원, 반례, 호환성, 과학적 주장, scope.
   검토 의견과 수용/거부 이유를 아래 기록하고 나서 구현한다.
3. C1–C5를 코드/metadata/회귀 테스트로 구현; 데이터 정규화나 모델 fit을
   기존 의미와 몰래 바꾸지 않고 schema/version 및 migration note를 기록한다.
4. 분석적으로 정답을 아는 synthetic + 전체 회귀 테스트; Windows/Linux와
   지원 Python 최소/현재 버전. 임계값을 완화해 통과시키지 않는다.
5. GitHub에서 기존 실제 데이터 1개를 selective LFS로 실행, corrected 필드와
   unknown/eligibility 경고가 산출물까지 보존되는지 검사한다. 이 데이터로
   threshold를 튜닝하거나 같은 데이터로 성능을 주장하지 않는다.
6. 통과한 PR만 병합하고 정확한 commit/CI 링크·잔여 실험 항목을 기록한다.

## 4. 코드로 완료할 수 없는 교정 실험 (미실행)

- 센서 calibration·noise floor, nominal capacity 정의, cycler counter 사양 확인.
- 동일 chemistry/SOC/protocol/time scale의 다중 온도 저항 측정으로 Ea의
  적합성·불확실성·온도 적용 범위를 추정; 별도 cell/batch hold-out 검증.
- matched endpoints 및 충분한 relaxation을 갖춘 전압 profile로 SOC alignment 확인.
- 독립 coulomb/energy 적분 reference와 추출값 비교, sampling gap 민감도 확인.
- full-cell proxy는 적절한 half-cell/OCV/EIS 등 독립 측정과 비교하되
  동일 데이터의 feature agreement를 독립적인 인과 증거로 세지 않는다.
- 오차 허용치·필요 표본 수·false-positive/abstention 목표는 실험 전에
  도메인 책임자가 정한다. 실제 결과 없이 calibrated로 표기하지 않는다.

외부 문헌 원문 접근은 현재 네트워크에서 차단됐다. 따라서 문헌 검증 완료나
새 DOI 근거를 주장하지 않는다. 위 식은 명시한 단위·정의에서 유도한 계약이며,
원문/장비 매뉴얼 대조 및 교정 데이터 승인은 별도 미완료 항목이다.

## 5. 독립 검토 기록

구현 전에 별도 검토자가 소스와 계획을 읽고 C1–C5를 조건부 승인했다.
아래 조건을 모두 수용한다. 아직 구현/실험 성공을 의미하지 않는다.

- C1: 중복 capacity의 상충 voltage를 임의 선택하지 않는다. acquisition 순서에서
  검증하고 band 경계를 적분 grid에 넣는다. 정규화 영역 적분 단위는 V, Wh가 아니다.
- C2: Ea는 양수만 허용하는 보편 상수가 아니라 명시적인 유한 signed fit 계수다.
  source/target 온도 모두 적용 범위 안에 있어야 하며 R(T)=A exp(Ea/RgT) 모델명을
  기록한다. calibration id, measurement definition, SOC/protocol/time scale 및
  온도 출처를 요구하되 실제 calibration 검증 완료는 주장하지 않는다.
- C3: total/CE를 legacy counter-derived로 보존하고 observed quantities를 별도 제공.
  observed VE는 독립적으로 (Ed/ΔQd)/(Ec/ΔQc)로 정의한다. 기존/새 feature의
  semantics 버전을 baseline·delta·diagnosis까지 전달하고 모르면 evidence에서 제외.
  중간 step을 건너뛰어 합쳐진 leg는 연속 관측인 것처럼 적분하지 않는다.
- C4: evidence/mode별 필요 측정 그룹·protocol·fit·baseline/feature semantics를
  명시한다. irrelevant rest/pulse 누락은 비관련 모드를 차단하지 않으며 필요한
  항목 누락은 높은 종합 score로 통과시키지 않는다. boolean truthiness를 제거하고
  per-mode eligibility/reason을 export에 보존한다. 새 confidence penalty는 만들지 않는다.
- C5: explicit A 또는 explicit Q_nom×C 기준을 공통 resolver로 전달하고 provenance
  기록. step/방향/clock reset/결측 경계에서 분리. 첫 연속 pulse만 count/stability에
  사용. 시간이 없으면 index를 초로 가정하지 않는다. 저항 proxy가 진짜 pre-pulse
  ΔV/ΔI 측정으로 검증되었다고 주장하지 않는다.

거부한 초기 제안(C-rate 용량 배율, CE>100 과충전 단정, chemistry 충돌 완화,
counter smoothing, 임의 quality 확률화)은 계속 제외한다.
검토의 acceptance tests(가역/상수-offset curve, signed Ea, observed offset 불변,
불연속 적분 거부, mode별 evidence 행렬, legacy/export roundtrip)를 회귀로 구현한다.