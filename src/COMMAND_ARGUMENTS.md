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

토큰은 `+` 또는 `,`로 조합할 수 있습니다. `-`는 내부적으로 `_`로 바뀌므로 `phase2-linear-solve`와 `phase2_linear_solve`는 같은 의미로 처리됩니다.

```text
default
phase1_topological_offset
phase2_linear_solve
phase2_linear_solve_collision_reject
phase2_newton_solve
phase2_qem_original
phase1_topological_offset+phase2_linear_solve
phase2_linear_solve+boundary_rail
config.json
```

## Supported Preset Tokens

| 토큰 | 적용 대상 | 효과 |
|---|---|---|
| `default` | 전체 | 기존 Phase 1 절차를 그대로 사용합니다. 모든 tetrahedron을 1→12 subdivision하고 non-adjacent tetrahedron을 제거하는 과정을 2회 수행합니다. |
| `phase1_topological_offset` | phase1 | Zint et al.의 simplicial embedding과 offset insertion을 수행한 뒤 기존 `retrieveCage()`로 boundary를 추출합니다. `topological_offset`도 같은 별칭으로 사용할 수 있습니다. |
| `boundary_rail` | phase1.5/phase2 | Build source-boundary anchors and closed cage edge loops, then constrain rail collapses to the ruled half-strips defined by source boundary tangents and outward co-normals. Must be combined with one of the four energy Phase 2 presets. |
| `phase2_linear_solve` | phase2 | `phase2Mode = "linear_solve"`; solve the QEM/quality linear system, accept valid points directly, and otherwise use Armijo backtracking from tangential smoothing |
| `phase2_linear_solve_collision_reject` | phase2 | Same defaults as `phase2_linear_solve` plus `phase2LinearSolveCollisionReject = true`; reject an invalid raw linear-solve point instead of backtracking |
| `phase2_newton_solve` | phase2 | `phase2Mode = "newton_solve"`; use Newton placement for every popped collapse candidate |
| `phase2_qem_original` | phase2 | `phase2Mode = "qem_original"`; pure Garland-Heckbert QEM edge-collapse cost/placement with non-QEM energies and collision rejection disabled |

`phase1_topological_offset`을 명시한 경우에만 새 Phase 1을 사용합니다. `default`와 기존 `phase2_*` 토큰만 사용한 명령은 이전과 동일하게 subdivision 기반 Phase 1을 수행합니다.

알 수 없는 토큰이 들어오면 `unknown parameter token` 오류와 함께 실행이 중단됩니다.

Newton placement 세부 가중치는 JSON에서 조절합니다. 주요 키는
`paramCageSimplifier.paramCollapse.phase2PlacementStrategy`,
`qemWeight`, `triangleQualityWeight`,
`uniformityWeight`, `newtonMaxIter`,
`lineSearchMaxIter`,
`phase2LinearSolveCollisionReject`입니다.

## JSON Config

JSON 파일을 첫 번째 인수로 넘기면 `ParamCageGenerator` 설정을 덮어씁니다. 예시는 `config.json`에 있습니다.

```json
{
  "paramCageInitializer": {
    "phase1Mode": "topological_offset"
  },
  "paramCageSimplifier": {
    "maxIter": 50,
    "relaxErrorIterStep": 5,
    "maxErrorRelaxIter": 4,
    "initError": 0.005,
    "errorStep": 0.005,
    "paramCollapse": {
      "maxValence": 8
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
| `paramCageInitializer.phase1Mode` | Phase 1 방식: `"subdivision"` 또는 `"topological_offset"` | `"subdivision"` |
| `paramCageSimplifier.maxIter` | cage simplifier 최대 반복 횟수 | `30` |
| `paramCageSimplifier.relaxErrorIterStep` | 에러 완화 반복 간격 | `5` |
| `paramCageSimplifier.maxErrorRelaxIter` | 최대 에러 완화 단계 | `4` |
| `paramCageSimplifier.initError` | 초기 Hausdorff distance 허용값 | `0.005` |
| `paramCageSimplifier.errorStep` | 에러 완화 단계별 증가값 | `0.005` |
| `paramCageSimplifier.paramCollapse.maxValence` | collapse 단계 최대 valence | `8` |
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

Simplicial embedding과 offset insertion으로 initial cage 생성:

```powershell
.\exeCageGenerator.exe phase1_topological_offset C:\models\bunny.obj C:\out 500
```

새 Phase 1과 Phase 2 preset 조합:

```powershell
.\exeCageGenerator.exe phase1_topological_offset+phase2_linear_solve C:\models\bunny.obj C:\out 500
```

nested cage 생성:

```powershell
.\exeCageGenerator.exe default C:\models\bunny.obj C:\out 1000 500 250
```

프리셋 조합 사용:

```powershell
.\exeCageGenerator.exe phase2_linear_solve+boundary_rail C:\models\bunny.obj C:\out 500
```

JSON 설정 파일 사용:

```powershell
.\exeCageGenerator.exe .\config.json C:\models\bunny.obj C:\out 500
```

## Output

출력은 기본적으로 `<output_dir_path>\<input_file_name_without_ext>\` 아래에 새 실행별 하위 폴더로 생성됩니다.

```text
<output_dir_path>\<input_name>\
  <run_timestamp>[__phase1_topological_offset]__phase2_<mode>[__<mode-specific-details>]\
```

기존 출력 경로와의 호환성을 위해 `default` 모드에는 Phase 1 접미사를 추가하지 않으며, 기존 `__collapse_hausdorff__flip_<mode>__relocate_<mode>` 접미사도 유지합니다. 새 모드에만 `__phase1_topological_offset`이 추가됩니다.

같은 하위 폴더명이 이미 있으면 `_001`, `_002`처럼 번호를 붙여 기존 결과를 덮어쓰지 않습니다. 날짜는 실행 시점의 시스템 날짜/시간이며 `YYYYMMDD_HHMMSS` 형식입니다.

| 파일 | 설명 |
|---|---|
| `log.txt` | 실행 로그입니다. |
| `<input_name>_debug_topological_offset.obj` | 새 Phase 1이 생성한 offset-inserted tetrahedral mesh의 face dump입니다. |
| `<input_name>_debug_retrieve_cage.obj` | Phase 1 boundary extraction 직후의 initial cage입니다. |
| `<input_name>_cage_0.obj` | 첫 번째 cage 결과입니다. |
| `<input_name>_cage_1.obj` | 두 번째 nested cage 결과입니다. 목표 정점 수를 여러 개 넘긴 경우 생성됩니다. |

## Notes

- `<output_dir_path>`는 실행 전에 존재해야 합니다. 프로그램은 그 아래에 입력 파일 이름의 하위 디렉터리와 실행별 결과 디렉터리를 만듭니다.
- `<target_Nv_i>`는 정수로 파싱됩니다. 의도와 다른 결과를 피하려면 양의 정수를 사용하세요.
- 입력 메쉬가 non-manifold이거나 non-watertight이면 경고를 남기고 계속 진행합니다. 단, 읽기 실패 또는 vertex/face가 없는 입력은 `invalid mesh`로 중단됩니다.
- 현재 `main.cpp`의 도움말 출력은 인덱스 표기가 일부 어긋나 있습니다. 실제 파싱 기준은 이 문서의 명령 형식입니다.
