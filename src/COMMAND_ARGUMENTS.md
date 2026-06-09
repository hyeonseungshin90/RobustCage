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

토큰은 `+` 또는 `,`로 조합할 수 있습니다. `-`는 내부적으로 `_`로 바뀌므로 `length-quality`와 `length_quality`는 같은 의미로 처리됩니다.

```text
default
collapse+flip_triangle_quality_hard+relocate_triangle_quality_hard
length+flip_triangle_quality_hard+relocate_triangle_quality_hard
length_quality_relative,flip_triangle_quality_hard
config.json
```

## Supported Preset Tokens

| 토큰 | 적용 대상 | 효과 |
|---|---|---|
| `default` | 전체 | 기본 파라미터를 그대로 사용합니다. |
| `collapse` / `collapse_default` / `collapse_hausdorff` | collapse | `priorityMode = "hausdorff"` |
| `length` | collapse | `priorityMode = "length"` |
| `length_quality` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "weighted"` |
| `length_quality_weighted` | collapse | `length_quality`와 동일합니다. |
| `length_quality_relative` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "relative_reject"` |
| `length_quality_absolute` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "absolute_reject"` |
| `length_quality_lexicographic` | collapse | `priorityMode = "length_quality"`, `lengthQualitySubMode = "lexicographic"` |
| `post_edge_length` | collapse | `priorityMode = "post_edge_length"` |
| `post_edge_length_hard` | collapse | `priorityMode = "post_edge_length_hard"` |
| `post_face_area` | collapse | `priorityMode = "post_face_area"` |
| `post_face_area_hard` | collapse | `priorityMode = "post_face_area_hard"` |
| `triangle_quality` | collapse | `priorityMode = "triangle_quality"` |
| `triangle_quality_hard` | collapse | `priorityMode = "triangle_quality_hard"` |
| `flip_valence` | flip | `priorityMode = "valence"` |
| `flip_triangle_quality_hard` | flip | `priorityMode = "triangle_quality_hard"` |
| `relocate_hausdorff` | relocate | `priorityMode = "hausdorff"` |
| `relocate_triangle_quality_hard` | relocate | `priorityMode = "triangle_quality_hard"` |

알 수 없는 토큰이 들어오면 `unknown parameter token` 오류와 함께 실행이 중단됩니다.

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
.\exeCageGenerator.exe length_quality_relative+flip_triangle_quality_hard C:\models\bunny.obj C:\out 500
```

JSON 설정 파일 사용:

```powershell
.\exeCageGenerator.exe .\config.json C:\models\bunny.obj C:\out 500
```

## Output

출력은 `<output_dir_path>\<input_file_name_without_ext>\` 아래에 생성됩니다.

| 파일 | 설명 |
|---|---|
| `log.txt` | 실행 로그입니다. |
| `<input_name>_cage_0.obj` | 첫 번째 cage 결과입니다. |
| `<input_name>_cage_1.obj` | 두 번째 nested cage 결과입니다. 목표 정점 수를 여러 개 넘긴 경우 생성됩니다. |

## Notes

- `<output_dir_path>`는 실행 전에 존재해야 합니다. 프로그램은 그 아래에 입력 파일 이름의 하위 디렉터리를 만듭니다.
- `<target_Nv_i>`는 정수로 파싱됩니다. 의도와 다른 결과를 피하려면 양의 정수를 사용하세요.
- 입력 메쉬가 non-manifold이거나 non-watertight이면 경고를 남기고 계속 진행합니다. 단, 읽기 실패 또는 vertex/face가 없는 입력은 `invalid mesh`로 중단됩니다.
- 현재 `main.cpp`의 도움말 출력은 인덱스 표기가 일부 어긋나 있습니다. 실제 파싱 기준은 이 문서의 명령 형식입니다.
