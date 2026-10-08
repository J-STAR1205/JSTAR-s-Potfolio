__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M03 2차 라운드 Step 5 (지침 25절) — 접합부 국소(junction-local) mesh 생성
#
# 분류: B. Parameter sensitivity / numerical convergence study.
#   dqdfdsoi.geo/.xao 원본은 전혀 수정하지 않는다 (gmsh.open()으로 읽기만 하고,
#   Field/generate 결과는 새 .msh 파일로만 저장).
#
# [NUMERICAL] [설계] QTCAD의 PoissonSolverParams.refined_region/h_refined은 기존에
#   정의된 Physical Volume 이름 단위로만 동작한다 — "접합부 ±5~10nm"만 떼어낸
#   전용 physical group이 원본 geometry에 없으므로(1차 라운드에서 확인), 그 전체
#   channel/source/drain 볼륨 전체를 h_refined 대상으로 쓸 수밖에 없었고, 이는
#   접합부보다 훨씬 넓은 영역(channel 길이 90nm)이라 h_refined=0.3에서 노드수가
#   2,673,597개까지 폭주했다(M03_mesh_convergence.py 1차 라운드, ArrayMemoryError).
#
#   이번에는 QTCAD의 adaptive refined_region 메커니즘을 전혀 쓰지 않고, gmsh의
#   네이티브 Mesh.Field(Box) 기능으로 "y=-45nm, y=+45nm 주변 ±half_width만" 국소
#   목표 크기(h_local)를 지정한다. Field는 기존 Physical Surface/Volume 번호를
#   전혀 바꾸지 않으므로 BooleanFragments로 생성된 원본 geometry의 번호 체계에
#   안전하다.
#
# [WARNING] 생성된 메시의 전체 노드 수를 QTCAD Poisson 파이프라인에 넣기 전에
#   먼저 확인한다 (지침 10절 + 1차 라운드 폭주 교훈 — 비용이 큰 단계 전에 먼저
#   값싼 진단을 한다).

import sys
import pathlib
import gmsh

script_dir = pathlib.Path(__file__).parent.resolve()
ascii_staging = pathlib.Path(
    r"C:\Users\norma\AppData\Local\Temp\claude\C--Users-norma-Desktop------QTCAD-Simulation"
    r"\1c089a01-c171-407b-9be8-42b4f3db9e1e\scratchpad\ascii_geo"
)
path_xao = ascii_staging / "dqdfdsoi.xao"  # [BASELINE] 기존 스크립트들과 동일한 원본(읽기 전용)

# [PHYSICS] 접합 위치: dqdfdsoi.geo의 치수 정의로부터 유도 (암기 아님)
#   channel_len = 6*5(gap) + 3*10(barrier) + 2*15(plunger) = 90nm -> channel: y in [-45,45]
#   즉 source/channel 접합은 y=-45nm, channel/drain 접합은 y=+45nm.
JUNCTION_LEFT_NM = -45.0
JUNCTION_RIGHT_NM = 45.0
HALF_WIDTH_NM = 10.0   # [PARAMETER] 접합 ±10nm만 국소 정제 (지침이 제시한 ±5~10nm 범위)
X_MARGIN_NM = 25.0     # domain_width=60 -> x in [-30,30]; 채널(폭 40, x in[-20,20])을
                        # 여유있게 덮도록 ±25nm 사용
Z_MIN_NM = -22.0        # buried_oxide 바닥까지 여유있게 포함 (box_thick+film_thick=20)
Z_MAX_NM = 3.0          # gate_oxide 상단까지 여유있게 포함 (gate_oxide_thick=2)
VOUT_NM = 100.0         # [NUMERICAL] box 바깥에서는 gmsh 기본/자동 크기가 그대로 쓰이도록
                        # 충분히 큰 값(장치 전체 크기보다 큼)을 지정 — Field는 다른 크기
                        # 제약과 min()으로 결합되므로 VOut이 크면 box 밖 영향 없음.
THICKNESS_NM = 5.0      # box 경계에서 h_local -> VOut으로 매끄럽게 전이


def generate_mesh(h_local_nm, out_msh_path):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.open(str(path_xao))

    fields = []
    for junction_y in (JUNCTION_LEFT_NM, JUNCTION_RIGHT_NM):
        f = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(f, "VIn", h_local_nm)
        gmsh.model.mesh.field.setNumber(f, "VOut", VOUT_NM)
        gmsh.model.mesh.field.setNumber(f, "XMin", -X_MARGIN_NM)
        gmsh.model.mesh.field.setNumber(f, "XMax", X_MARGIN_NM)
        gmsh.model.mesh.field.setNumber(f, "YMin", junction_y - HALF_WIDTH_NM)
        gmsh.model.mesh.field.setNumber(f, "YMax", junction_y + HALF_WIDTH_NM)
        gmsh.model.mesh.field.setNumber(f, "ZMin", Z_MIN_NM)
        gmsh.model.mesh.field.setNumber(f, "ZMax", Z_MAX_NM)
        gmsh.model.mesh.field.setNumber(f, "Thickness", THICKNESS_NM)
        fields.append(f)

    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)

    gmsh.model.mesh.generate(3)
    gmsh.write(str(out_msh_path))

    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    n_elem = sum(len(tags) for dim, tag in gmsh.model.mesh.getElements()[1]
                 for tags in [tag] if dim == 3) if False else None
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=3)
    n_tet = sum(len(t) for t in elem_tags)
    gmsh.finalize()
    return n_nodes, n_tet


if __name__ == "__main__":
    h_values_nm = [2.0, 1.0, 0.5]
    print(f"{'h_local(nm)':>12} {'nodes':>12} {'tets':>12}")
    for h in h_values_nm:
        out_path = ascii_staging / f"dqdfdsoi_junction_local_h{h}.msh"
        n_nodes, n_tet = generate_mesh(h, out_path)
        print(f"{h:>12.2f} {n_nodes:>12d} {n_tet:>12d}  -> {out_path.name}")
