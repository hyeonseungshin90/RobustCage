# Command Arguments

`exeCageGenerator.exe`는 위치 기반 명령 인수를 사용합니다. 최소 1개의 목표 정점 수가 필요합니다.

```text
exeCageGenerator.exe <parameters_or_config> <input_model_path> <output_dir_path> <target_Nv_0> [target_Nv_1 ... target_Nv_n]
```

## Arguments

| 순서 | 이름 | 필수 | 설명 |
|---:|---|:---:|---|
| 0 | `exeCageGenerator.exe` | 예 | 실행 파일 이름입니다. CMake 대상 이름은 `exeCageGenerator`입니다. |
| 1 | `<parameters_or_config>` | 예 | `default`, 프리셋 토큰, 토큰 조합, 또는 JSON 설정 파일 경로입니다. |
| 2 | `<input_model_path>` | 예 | 입력 메쉬 파일 경로입니다. 코드에서 일반 파일인지 검사합니다. |
| 3 | `<output_dir_path>` | 예 | 출력 디렉터리 경로입니다. 이미 존재하는 디렉터리여야 합니다. |
| 4+ | `<target_Nv_i>` | 예 | 생성할 cage의 목표 정점 수입니다. 여러 개를 넘기면 nested cage를 순서대로 생성합니다. |

## Parameter Argument

첫 번째 인수 `<parameters_or_config>`는 두 가지 방식으로 해석됩니다.

1. 같은 경로에 실제 파일이 있으면 JSON 설정 파일로 읽습니다.
2. 파일이 아니면 프리셋 토큰 문자열로 처리합니다.

토큰은 `+` 또는 `,`로 조합할 수 있습니다. `-`는 내부적으로 `_`로 바뀌므로 `collapse-length-quality`와 `collapse_length_quality`는 같은 의미로 처리됩니다.

```text
default
collapse+flip_triangle_quality_hard+relocate_triangle_quality_hard
collapse_length+flip_triangle_quality_hard+relocate_triangle_quality_hard
collapse_length_quality_relative,flip_triangle_quality_hard
config.json
```

## Supported Preset Tokens

| 토큰 | 적용 대상 | 효과 |
|---|---|---|
| `default` | 전체 | 기본 파라미터를 그대로 사용합니다. |
| `collapse` / `collapse_default` / `collapse_hausdorff` | collapse | `priorityMode = "hausdorff"` |
| `collapse_length` | collapse | `priorityMode = "length"` |
| `collapse_length_quality` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "weighted"` |
| `collapse_length_quality_weighted` | collapse | `collapse_length_quality`와 동일합니다. |
| `collapse_length_quality_relative` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "relative_reject"` |
| `collapse_length_quality_absolute` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "absolute_reject"` |
| `collapse_length_quality_lexicographic` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "lexicographic"` |
| `collapse_post_edge_length` | collapse | `priorityMode = "post_edge_length"` |
| `collapse_post_edge_length_hard` | collapse | `priorityMode = "post_edge_length_hard"` |
| `collapse_post_face_area` | collapse | `priorityMode = "post_face_area"` |
| `collapse_post_face_area_hard` | collapse | `priorityMode = "post_face_area_hard"` |
| `collapse_triangle_quality` | collapse | `priorityMode = "triangle_quality"` |
| `collapse_triangle_quality_hard` | collapse | `priorityMode = "triangle_quality_hard"` |
| `phase2_newton` | phase2 | `phase2Mode = "newton"`, `collapsePlacementMethod = "optimization"`, `phase2PlacementStrategy = "adaptive"`, `curvatureMode = "weighted_qem"`, `uniformityMode = "source"`, `positionFidelityWeight = 1.0`, barrier weights `0.0` by default because exact checks are hard constraints |
| `phase2_adaptive` | phase2 | QEM linear solve first, Newton only when local quality/residual/final-collapse criteria request it |
| `phase2_linear_only` | phase2 | `phase2Mode = "linear_only"`; QEM linear solve placement during Phase 2 collapse with no Newton refinement |
| `phase2_qem` | phase2 | `phase2Mode = "qem"`; for a watertight manifold cage, estimate the target mean edge length as `initialMean * sqrt((initialVertices - chi) / (targetVertices - chi))`, where `chi = V - E + F`. Queue ranking combines QEM with `current edge length / target mean edge length`, so shorter edges receive higher priority. The previous source-adaptive queue ranking remains available by combining `phase2_qem+uniformity_source`. Accept valid QEM points directly, otherwise backtrack from tangential smoothing; skip the collapse candidate if no backtracking point satisfies hard intersection/validity checks |
| `phase2_qem_reject_on_collision` / `qem_reject_on_collision` | phase2 | `phase2Mode = "qem"` (same defaults as `phase2_qem`) plus `phase2QemRejectOnLinearCollision = true`. If the raw QEM linear-solve point fails hard validity (collision/degenerate/wrinkle), reject the edge outright instead of backtracking from tangential smoothing. Standalone preset — unlike `phase2_qem`, no need to combine tokens. `phase2_qem` alone keeps the default `false` (backtracking) behavior. |
| `phase2_qem_no_collision` | phase2 | `phase2Mode = "qem_no_collision"`; pure Garland-Heckbert QEM edge-collapse cost/placement only, with non-QEM energies and collision rejection disabled |
| `phase2_final_newton` | phase2 | QEM linear solve placement during collapse, then one fixed-topology Newton relocation pass |
| `phase2_newton_only` | phase2 | Skip QEM linear solve; use Newton placement for every popped collapse candidate |
| `phase2_quadratic_surrogate` | phase2 | Separate experimental strategy: keep QEM queue ranking, then add the 4x4 quadratic surrogate as an extra placement candidate with QEM/fallback safeguards; no Newton refinement |
| `flip_valence` | flip | `priorityMode = "valence"` |
| `flip_triangle_quality_hard` | flip | `priorityMode = "triangle_quality_hard"` |
| `relocate_hausdorff` | relocate | `priorityMode = "hausdorff"` |
| `relocate_triangle_quality_hard` | relocate | `priorityMode = "triangle_quality_hard"` |
| `collapse_optimization` / `collapse_newton` | phase2/collapse | `phase2Mode = "newton"`, `collapsePlacementMethod = "optimization"`, `curvatureMode = "weighted_qem"`, `uniformityMode = "source"` |
| `collapse_sampling` | collapse | `collapsePlacementMethod = "sampling"` |
| `newton_damped` | collapse | `newtonSolverMode = "damped"` |
| `newton_trust_region` | collapse | `newtonSolverMode = "trust_region"` |
| `robust_exact_reject` | collapse | `robustnessMode = "exact_reject"` |
| `robust_exact_backtracking` | collapse | `robustnessMode = "exact_backtracking"` |
| `robust_ipc` | collapse | `robustnessMode = "ipc_line_search"` |
| `curvature_none` | collapse | `curvatureMode = "none"` |
| `curvature_weighted_qem` | collapse | `curvatureMode = "weighted_qem"` |
| `curvature_normal_matching` | collapse | `curvatureMode = "normal_matching"` |
| `uniformity_none` | collapse | `uniformityMode = "none"` |
| `uniformity_source` | collapse | `uniformityMode = "source"`; preserve the source-adaptive target-edge calculation |
| `uniformity_global` | collapse | `uniformityMode = "global"`; rank by current cage edge length relative to the Euler-estimated target mean edge length |

알 수 없는 토큰이 들어오면 `unknown parameter token` 오류와 함께 실행이 중단됩니다.

Newton placement 세부 가중치는 JSON에서 조절합니다. 주요 키는
`paramCageSimplifier.paramCollapse.phase2PlacementStrategy`,
`qemWeight`, `selfBarrierWeight`,
`originalBarrierWeight`, `positionFidelityWeight`, `curvatureWeight`, `triangleQualityWeight`,
`uniformityWeight`, `barrierActivationDistanceFactor`, `newtonMaxIter`,
`lineSearchMaxIter`, `lineSearchCcdSamples`,
`phase2NewtonQualityThreshold`, `phase2NewtonResidualThreshold`,
`phase2NewtonResidualGrowth`, `phase2NewtonFinalRefineCollapses`,
`phase2QemRejectOnLinearCollision`입니다.

## JSON Config

JSON 파일을 첫 번째 인수로 넘기면 `ParamCageGenerator` 설정을 덮어씁니다. 예시는 `config.json`에 있습니다.

```json
{
  "paramCageSimplifier": {
    "maxIter": 50,
    "relaxErrorIterStep": 5,
    "maxErrorRelaxIter": 4,
    "initError": 0.005,
    "errorStep": 0.005,
    "paramCollapse": {
      "maxValence": 8,
      "priorityMode": "hausdorff"
    },
    "paramFlip": {
      "maxValence": 8
    }
  }
}
```

주요 JSON 키는 다음과 같습니다.

| 경로 | 설명 | 기본값 |
|---|---|---:|
| `paramCageSimplifier.maxIter` | cage simplifier 최대 반복 횟수 | `30` |
| `paramCageSimplifier.relaxErrorIterStep` | 에러 완화 반복 간격 | `5` |
| `paramCageSimplifier.maxErrorRelaxIter` | 최대 에러 완화 단계 | `4` |
| `paramCageSimplifier.initError` | 초기 Hausdorff distance 허용값 | `0.005` |
| `paramCageSimplifier.errorStep` | 에러 완화 단계별 증가값 | `0.005` |
| `paramCageSimplifier.paramCollapse.maxValence` | collapse 단계 최대 valence | `8` |
| `paramCageSimplifier.paramCollapse.priorityMode` | collapse 우선순위 모드 | `"hausdorff"` |
| `paramCageSimplifier.paramCollapse.lengthQualitySubMode` | length-quality 세부 모드 | `"weighted"` |
| `paramCageSimplifier.paramCollapse.lengthQualityWeight` | length-quality 가중치 | `5.0` |
| `paramCageSimplifier.paramCollapse.lengthQualityDegradationRatio` | relative reject 품질 저하 허용 비율 | `0.5` |
| `paramCageSimplifier.paramCollapse.lengthQualityMinQuality` | absolute reject 최소 triangle quality | `0.1` |
| `paramCageSimplifier.paramRelocate.smoothIter` | relocate smoothing 반복 횟수 | `3` |
| `paramCageSimplifier.paramRelocate.priorityMode` | relocate 우선순위 모드 | `"hausdorff"` |
| `paramCageSimplifier.paramFlip.maxValence` | flip 단계 최대 valence | `8` |
| `paramCageSimplifier.paramFlip.priorityMode` | flip 우선순위 모드 | `"valence"` |

주의: JSON 역직렬화에서 일부 필드는 필수로 접근합니다. `paramCageSimplifier`, `maxIter`, `relaxErrorIterStep`, `maxErrorRelaxIter`, `initError`, `errorStep`, `paramCollapse`, `paramCollapse.maxValence`, `paramFlip.maxValence`는 예시처럼 포함하는 것이 안전합니다.

## Examples

단일 cage 생성:

```powershell
.\exeCageGenerator.exe default C:\models\bunny.obj C:\out 500
```

nested cage 생성:

```powershell
.\exeCageGenerator.exe default C:\models\bunny.obj C:\out 1000 500 250
```

프리셋 조합 사용:

```powershell
.\exeCageGenerator.exe collapse_length_quality_relative+flip_triangle_quality_hard C:\models\bunny.obj C:\out 500
```

JSON 설정 파일 사용:

```powershell
.\exeCageGenerator.exe .\config.json C:\models\bunny.obj C:\out 500
```

## Output

출력은 기본적으로 `<output_dir_path>\<input_file_name_without_ext>\` 아래에 새 실행별 하위 폴더로 생성됩니다.

```text
<output_dir_path>\<input_name>\
  <run_timestamp>__collapse_<mode>__flip_<mode>__relocate_<mode>\
```

같은 하위 폴더명이 이미 있으면 `_001`, `_002`처럼 번호를 붙여 기존 결과를 덮어쓰지 않습니다. 날짜는 실행 시점의 시스템 날짜/시간이며 `YYYYMMDD_HHMMSS` 형식입니다.

| 파일 | 설명 |
|---|---|
| `log.txt` | 실행 로그입니다. |
| `<input_name>_cage_0.obj` | 첫 번째 cage 결과입니다. |
| `<input_name>_cage_1.obj` | 두 번째 nested cage 결과입니다. 목표 정점 수를 여러 개 넘긴 경우 생성됩니다. |

## Notes

- `<output_dir_path>`는 실행 전에 존재해야 합니다. 프로그램은 그 아래에 입력 파일 이름의 하위 디렉터리와 실행별 결과 디렉터리를 만듭니다.
- `<target_Nv_i>`는 정수로 파싱됩니다. 의도와 다른 결과를 피하려면 양의 정수를 사용하세요.
- 입력 메쉬가 non-manifold이거나 non-watertight이면 경고를 남기고 계속 진행합니다. 단, 읽기 실패 또는 vertex/face가 없는 입력은 `invalid mesh`로 중단됩니다.
- 현재 `main.cpp`의 도움말 출력은 인덱스 표기가 일부 어긋나 있습니다. 실제 파싱 기준은 이 문서의 명령 형식입니다.
