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
phase1_topological_offset+boundary_rail_compare
config.json
```

## Supported Preset Tokens

| 토큰 | 적용 대상 | 효과 |
|---|---|---|
| `default` | 전체 | 기존 Phase 1 절차를 그대로 사용합니다. 모든 tetrahedron을 1→12 subdivision하고 non-adjacent tetrahedron을 제거하는 과정을 2회 수행합니다. |
| `phase1_topological_offset` | phase1 | Zint et al.의 simplicial embedding과 offset insertion을 수행한 뒤 기존 `retrieveCage()`로 boundary를 추출합니다. `topological_offset`도 같은 별칭으로 사용할 수 있습니다. |
| `boundary_rail` | phase1.5/phase2 | Source-boundary anchor를 삽입하고 Dijkstra로 초기 loop를 구성한 뒤, anchor를 고정한 intrinsic flip geodesics와 실제 mesh split으로 rail 경로를 단축합니다. 이후 rail collapse는 source boundary tangent와 outward co-normal의 ruled half-strip으로 제약합니다. 네 가지 energy Phase 2 preset 중 하나와 조합해야 합니다. |
| `boundary_rail_vertex` | phase1.5/phase2 | `boundary_rail`과 동일하지만 ray를 boundary vertex에서 인접한 두 edge co-normal의 bisector 방향으로 쏩니다. 비교 실험용이며 source-edge support에는 평균 전의 edge co-normal을 저장합니다. Energy Phase 2 preset과 조합해야 합니다. |
| `boundary_rail_compare` | phase1/benchmark | Phase 1을 한 번만 실행한 뒤 동일 initial cage와 source의 독립 복사본에서 edge-midpoint 방식과 vertex-bisector 방식을 각각 실행합니다. simplification은 생략하고 run-local CSV와 두 결과 cage/rail을 기록합니다. |
| `phase2_linear_solve` | phase2 | `phase2Mode = "linear_solve"`; repeat linear-solve collapses with Armijo backtracking and quality-priority flips, then run final relocation sweeps combining Voronoi tangential smoothing with source-surface attraction |
| `phase2_linear_solve_collision_reject` | phase2 | Same defaults as `phase2_linear_solve` plus `phase2LinearSolveCollisionReject = true`; reject an invalid raw linear-solve point instead of backtracking |
| `phase2_newton_solve` | phase2 | `phase2Mode = "newton_solve"`; use Newton placement for every popped collapse candidate |
| `phase2_qem_original` | phase2 | `phase2Mode = "qem_original"`; pure Garland-Heckbert QEM edge-collapse cost/placement with non-QEM energies and collision rejection disabled |

`phase1_topological_offset`을 명시한 경우에만 새 Phase 1을 사용합니다. `default`와 기존 `phase2_*` 토큰만 사용한 명령은 이전과 동일하게 subdivision 기반 Phase 1을 수행합니다.

`linear_solve`에서는 `collapse → flip → rail update`를 기본 최대 30회 반복한 뒤,
마지막에 relocation 단계를 한 번 실행합니다.
Flip과 rail 라벨 갱신 이후 collapse 후보 큐를 다시 만들어 새롭게 가능한 collapse를
재시도하며, global target edge length는 최초 값으로 유지합니다. 목표 정점 수에
도달하면 collapse를 생략하고 flip과 rail update는 계속할 수 있습니다.
Rail update는 boundary rail이 활성화된 경우에만 실행합니다. Collapse, flip,
rail update가 모두 없으면 반복을 조기 종료합니다. 반복 종료 후 relocation을 수행하며,
그 뒤에는 collapse, flip, rail update를 다시 실행하지 않습니다. Flip은 boundary rail
활성화 여부와 관계없이 두 삼각형의 최소 quality 증가량이 큰 순서로 수행합니다.
교차가 생기는 후보는 reject합니다.
Collapse 위치의 linear solve와 line-search energy는 plane과 triangle-quality
surrogate를 사용하며 uniformity는 제외합니다. Uniformity는 collapse 우선순위에만
적용하며, `global` 모드에서는 `uniformityWeight * (edge_length / target_length)^2`를
큐 점수에 더합니다.
Relocation은 이웃 정점의 mixed Voronoi 면적을 가중치로 사용하는
Botsch–Kobbelt tangential smoothing point `t`와 현재 정점에서 원본 표면의
최근접점 `s`를 결합합니다. 두 목표점은 해당 정점의 line search 동안 고정하며,
`E(x) = tangentialWeight * ||x-t||² + surfaceWeight * ||x-s||²`를 사용합니다.
현재 위치에서 두 목표점의 가중 평균 방향으로
`alpha = 1, 1/2, 1/4, ...` backtracking하고, 이 에너지가 엄격히 감소하며
비퇴화·비반전 및 교차 검사를 모두 통과해야 수락합니다. 주변 삼각형의 최소
quality는 작은 수치 오차 허용 범위 내에서
`new_min_quality >= min(old_min_quality, minTriangleQuality)`를 만족해야 합니다.
Quality 범위는 `[0, 1]`이고 기본 하한은 `0.2`입니다. 기존 quality가 하한보다
높으면 하한까지 낮아질 수 있지만, 이미 하한보다 낮으면 더 나빠질 수 없습니다.
마지막 relocation 단계는 전체 최대 20 sweep을 수행하며, 매 sweep에서 현재 위치를 기준으로
Voronoi 면적·normal·원본 최근접점을 재계산합니다. 한 sweep에서 이동이 없으면 조기
종료합니다. 정점별 line search는 각 sweep에서 최대 12회 시도합니다.
`+boundary_rail`에서 한 끝점만 rail에 속하는 mixed edge는 일반 정점을 기존 rail
정점으로 합치는 방향으로만 collapse할 수 있습니다. 새 위치는 기존 rail 정점의
위치로 고정하고 rail ID와 loop 연결을 유지하며, 기존 collapse 유효성 검사를
통과해야 합니다. 실제 rail edge의 collapse는 계속 source-boundary support
half-strip에 새 위치를 투영합니다.
Rail edge의 flip은 금지되며, rail 정점과 실제 mesh boundary 정점은 relocation
중 고정됩니다. Flip과 relocation은 정점 수를 유지합니다.

Linear-solve flip은 triangle quality로 우선순위와 개선 여부를 판단합니다.
모든 flip은 두 삼각형의 유한한 최소 quality 증가량이 `1e-12`보다 클 때만 허용합니다.
Rail chord 개수는 우선순위나 수락 여부에 영향을 주지 않습니다.
모든 후보에 기존 위상·valence·퇴화·source/cage 교차 검사를 적용하며,
별도의 flip quality 하한은 추가하지 않습니다. 갱신된 주변 후보는 다시 큐에 넣고
수락 직전에 재검사합니다. Flip은 rail 라벨·정점 위치·edge·둘레를 보존하지만
cage 표면은 바뀔 수 있습니다.
실제 rail edge collapse와 기존 projection 및 line search는 유지합니다.

C++ 진입점은 `do_quality_flip()`이며, chord 정책을 선택하는 인수나 별도의
chord-aware 모드는 없습니다. Boundary rail을 활성화하면 flip 뒤의 삼각형 단위
rail update를 실행합니다.

초기 rail 구성은 Phase 2 이전에 `Dijkstra → intrinsic flip geodesics → 실제 mesh split`
순서로 수행합니다. Source-boundary anchor를 cage에 삽입한 뒤 기존 Dijkstra 경로
탐색과 검증으로 단일 폐곡선들을 구성합니다. 이어서 Geometry Central의 flip geodesics로
전체 rail network를 단축하며 모든 anchor를 고정합니다. 단축된 경로를 원래 cage 표면에
trace하고, edge 교차점에 실제 정점을 삽입하고 주변 edge와 face를 분할하여 경로를
실제 rail edge들의 연결로 만듭니다. 이때 기존 정점 위치와 cage 표면은 유지되지만
정점과 삼각형 수는 증가할 수 있습니다. 접힌 면을 가로지르는 경로를 양 끝점 사이의
직선 3D edge로 대체하지 않습니다. 이후에는 기존 rail edge collapse를 적용합니다.

이 과정은 새 CLI 옵션 없이 `boundary_rail`과 `boundary_rail_vertex` 모두에 적용하며,
`boundary_rail_compare`의 두 독립 구성에도 적용합니다. Intrinsic 단축이나 실제 mesh
반영에 실패하면 Dijkstra 구성을 마친 시점의 cage와 rail을 유지하고 warning과 실패
통계를 기록합니다. 이 경우 결과에는 geodesic 갱신이 적용되지 않은 것입니다.
Flip geodesics는 locally shortest path를 구하며 전역 최단 경로를 보장하지 않습니다.
고정된 anchor 때문에 경로가 바뀌지 않을 수도 있습니다. 초기 rail 구성에서만 실행하며,
Phase 2 중 rail을 다시 geodesic으로 단축하지 않습니다. 이후 collapse가 초기 경로의
geodesic 성질을 유지한다는 보장도 없습니다.

Boundary rail을 활성화한 `linear_solve`에서는 매 quality-flip 단계 뒤에
삼각형 단위 `update_boundary_rails()`를 호출합니다. 이 함수는 더 이상 가능한
갱신이 없을 때까지 삼각형의 rail 경로 `A-B-C`를 `A-C`로 반복 교체하고 `B`의 라벨을
해제하며, 최소 3정점의 단일 폐곡선을 보존합니다. 라벨만 변경하는 함수이고
source 거리나 opening 크기에 대한 별도 오차 제한은 없습니다.
Rail update만 발생한 cycle도 진행으로 간주하므로, 다음 cycle에서 갱신된 제약으로
collapse와 flip을 다시 시도합니다. `phase2QualityPolishIterations = 0`은 기존
collapse-only 동작을 유지하여 flip, rail update, 최종 relocation을 생략합니다.

이 에너지는 전체 Hausdorff distance를 최소화하거나 그 개선을 보장하지 않습니다.
교차 검사는 후보 위치 기준이며, 이동 경로 전체의 CCD나 양의 표면 간격을
보장하는 검사는 아닙니다.
이 품질 개선 단계는 기존 `paramFlip.priorityMode`, `requireRegularValence`,
`paramRelocate.priorityMode`, `smoothIter` 대신 위의 energy·quality 기준과
`phase2QualityPolishIterations`, `paramRelocate.qualitySweeps`,
`paramRelocate.lineSearchMaxIter`, `paramRelocate.tangentialWeight`,
`paramRelocate.surfaceWeight`, `paramRelocate.minTriangleQuality`를 사용합니다.
두 가중치는 유한한 비음수여야 합니다. `surfaceWeight = 0`이면 tangential target만,
`tangentialWeight = 0`이면 원본 최근접점 target만 사용하며, 둘 다 `0`이면 relocation을
생략합니다. `minTriangleQuality`는 유한한 `[0, 1]` 값이어야 합니다. 기존 JSON 키 이름
`phase2QualityPolishIterations`는 호환성을 위해 유지하며, collapse/flip/rail-update 반복 한도를
뜻합니다. `0`이면 collapse만 한 번 실행하고 flip, rail update, 마지막 relocation을 모두 생략합니다.

알 수 없는 토큰이 들어오면 `unknown parameter token` 오류와 함께 실행이 중단됩니다.

Newton placement 세부 가중치는 JSON에서 조절합니다. 주요 키는
`paramCageSimplifier.paramCollapse.phase2PlacementStrategy`,
`planeWeight`, `triangleQualityWeight`,
`uniformityWeight`, `newtonMaxIter`,
`lineSearchMaxIter`,
`phase2LinearSolveCollisionReject`입니다.

`planeWeight` controls the local point-to-plane approximation energy. Legacy JSON
`qemWeight` remains accepted; `planeWeight` takes precedence when both keys are present.
`curvatureMode` uses `none` or `weighted_plane`. JSON values `weighted-plane`,
`weighted_qem`, and `weighted-qem` are accepted and normalized to `weighted_plane`.

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
| `paramCageSimplifier.maxIter` | 기존 일반 cage simplifier 최대 반복 횟수 | `30` |
| `paramCageSimplifier.phase2QualityPolishIterations` | linear-solve의 collapse/flip/rail-update 반복 한도. `0`이면 collapse만 한 번 실행하고 flip, rail update, 마지막 relocation 생략 | `30` |
| `paramCageSimplifier.relaxErrorIterStep` | 에러 완화 반복 간격 | `5` |
| `paramCageSimplifier.maxErrorRelaxIter` | 최대 에러 완화 단계 | `4` |
| `paramCageSimplifier.initError` | 초기 Hausdorff distance 허용값 | `0.005` |
| `paramCageSimplifier.errorStep` | 에러 완화 단계별 증가값 | `0.005` |
| `paramCageSimplifier.paramCollapse.maxValence` | collapse 단계 최대 valence | `8` |
| `paramCageSimplifier.paramRelocate.smoothIter` | 기존 일반 relocate smoothing 반복 횟수 | `3` |
| `paramCageSimplifier.paramRelocate.qualitySweeps` | linear-solve의 마지막 relocation 단계 전체 sweep 한도. `0`이면 relocation 생략 | `20` |
| `paramCageSimplifier.paramRelocate.lineSearchMaxIter` | linear-solve relocation energy의 정점당 backtracking 시도 한도 | `12` |
| `paramCageSimplifier.paramRelocate.tangentialWeight` | tangential smoothing target의 제곱거리 가중치. 유한한 비음수 | `1.0` |
| `paramCageSimplifier.paramRelocate.surfaceWeight` | 원본 최근접점 target의 제곱거리 가중치. 유한한 비음수. 두 가중치 모두 `0`이면 relocation 생략 | `1.0` |
| `paramCageSimplifier.paramRelocate.minTriangleQuality` | `[0, 1]`의 quality 하한. 기존 최소 quality가 이미 하한보다 낮으면 악화를 금지 | `0.2` |
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

동일한 Phase 1 cage에서 두 anchor 방법 비교(`target_Nv`는 CLI 형식상 필요하지만 사용하지 않음):

```powershell
.\exeCageGenerator.exe phase1_topological_offset+boundary_rail_compare C:\models\open.obj C:\out 0
```

이 모드는 `<run>/boundary_rail_anchor_benchmark.csv`에 detected, ray-valid,
built loop 수와 model별 winner를 쓰고, `<model>_initial_cage.obj`,
`<model>_edge_anchor_cage.obj`, `<model>_vertex_anchor_cage.obj` 및 두 rail 파일을 저장합니다.

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
| `<input_name>_cage_<i>_rails.obj` | cage `<i>`의 boundary rail만 담은 파일입니다. rail edge는 OBJ line element(`l`), rail vertex는 point element(`p`)로 기록되고 rail id마다 `o rail_<id>` / `g rail_<id>` 그룹으로 나뉩니다. rail이 하나도 없으면 생성되지 않습니다. |
| `<input_name>_cage_<i>_rails.txt` | 같은 rail을 cage OBJ의 1-based vertex index로 적은 목록입니다. `V <rail_id> <cage_vertex> <x> <y> <z>`와 `E <rail_id> <cage_vertex_a> <cage_vertex_b>` 두 종류의 줄로 구성되며, cage를 불러온 뒤 rail 정점/간선을 그 안에서 찾아 강조할 때 사용합니다. |

## Notes

- `<output_dir_path>`는 실행 전에 존재해야 합니다. 프로그램은 그 아래에 입력 파일 이름의 하위 디렉터리와 실행별 결과 디렉터리를 만듭니다.
- `<target_Nv_i>`는 정수로 파싱됩니다. 의도와 다른 결과를 피하려면 양의 정수를 사용하세요.
- rail 파일(`*_rails.obj`, `*_rails.txt`)은 `boundary_rail` 프리셋으로 실행해 rail이 실제로 만들어졌을 때만 생성됩니다. rail은 입력 메쉬의 boundary loop에서 나오므로 watertight 입력에서는 생성되지 않습니다.
- 입력 메쉬가 non-manifold이거나 non-watertight이면 경고를 남기고 계속 진행합니다. 단, 읽기 실패 또는 vertex/face가 없는 입력은 `invalid mesh`로 중단됩니다.
- 현재 `main.cpp`의 도움말 출력은 인덱스 표기가 일부 어긋나 있습니다. 실제 파싱 기준은 이 문서의 명령 형식입니다.
