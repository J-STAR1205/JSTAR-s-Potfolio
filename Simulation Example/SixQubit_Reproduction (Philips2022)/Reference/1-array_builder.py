__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [EXTENSION] Six-qubit Si/SiGe 선형 배열 재현 — Phase 0/1 (지오메트리 빌드)
#
# 레퍼런스: Philips et al., "Universal control of a six-qubit quantum processor
#   in silicon", Nature 609, 919 (2022). (PDF: 상위 "QTCAD Simulation" 폴더)
#
# 분류: C. 물리 확장(physics extension) — 공식 튜토리얼이 아님. 새 geometry.
#
# [PHYSICS] 논문에서 직접 확인한 수치 (Methods + Supplementary + Fig.1a SEM 이미지):
#   - interdot pitch = 90 nm (Device C, 본문 데이터)
#   - 게이트 토폴로지: SD1 - B0 P1 B1 P2 B2 P3 B3 P4 B4 P5 B5 P6 B6 - SD2
#   - 개별 plunger/barrier 폭은 본문에 숫자로 없음(figure 스케일바로만 추정 가능) —
#     여기서는 M03/M17 튜토리얼의 plunger:barrier 비율(15:10, 즉 3:2)을 유지하면서
#     합 + gap*2 = 90nm pitch에 맞게 선형 스케일: barrier_w=25, plunger_w=45, gap=10
#     (25+45+2*10=90). *** 이 치수는 1차 추정치이며, 실제 SEM 스케일바 측정으로
#     나중에 보정이 필요함을 명시 (설계도안 7절 리스크 항목) ***
#
# [NUMERICAL] N_dots를 파라미터화해서 Phase 1(N_dots=2, 기존 M03/M06 규모와 동등)로
#   먼저 파이프라인을 검증하고, Phase 2에서 동일 코드로 N_dots=6을 생성한다 — 전체를
#   한 번에 만들지 않는다(지침 원칙, M03 1차 라운드 폭주 교훈).
#
# [WARNING] sensing dot(SD1/SD2)과 EDSR 전용 screening gate, micromagnet은 이번
#   Phase에 포함하지 않는다 — 순수 plunger/barrier 체인의 전기정전학적 형성만
#   먼저 검증한다 (설계도안 5절 Phase 1-2).

import shutil
import sys
from pathlib import Path
from qtcad.builder import Builder, Polygon, Mask

N_DOTS = int(sys.argv[1]) if len(sys.argv) > 1 else 2  # [PARAMETER] Phase 1=2, Phase 2=6

script_dir = Path(__file__).parent.resolve()

# Length scales -----------------------------------------------------------------
char_len = 4  # characteristic mesh length (M03/M17과 동일 기본값)

domain_w = 60       # x폭 (M03/M17과 동일)
channel_w = 40      # 채널(실리콘) 폭
plunger_w = 45      # [PHYSICS, 1차 추정] plunger gate 길이 (y방향)
barrier_w = 25      # [PHYSICS, 1차 추정] barrier gate 길이
pitch = 10          # 게이트 간 gap
source_drain_w = 20  # source/drain(=SD1/SD2 placeholder) 길이, M03과 동일값 재사용

QD_w = plunger_w + barrier_w + 2 * pitch  # = 90nm, 논문의 interdot pitch와 일치하도록 설계
assert QD_w == 90, f"QD_w={QD_w}, expected 90nm interdot pitch"

# [NUMERICAL] 버그 수정: N_dots개 plunger 중심 간 거리(interdot pitch)는 QD_w로
# 정확하지만(plunger_w/2+gap+barrier_w+gap+plunger_w/2=QD_w), 체인 양 끝에는
# barrier가 "추가로 하나 더" 필요하다(N_dots개 plunger -> N_dots+1개 barrier).
# 따라서 총 길이는 N_dots*QD_w가 아니라 N_dots*QD_w + barrier_w (끝 barrier 1개분).
channel_len = N_DOTS * QD_w + barrier_w
domain_l = 2 * source_drain_w + channel_len

box_thick = 10
channel_thick = 10  # [PHYSICS] 논문의 QW+spacer+cap 스택(39nm)을 1차 근사로 단순화 —
                     # Phase 3(valley splitting)부터 실제 층 구조(cap 1nm/spacer 30nm/
                     # QW 8nm)로 세분화 예정. Phase 1-2는 M03/M17과 동일한 단층 근사.
EOT_thick = 2

print(f"[Phase {'1 (검증)' if N_DOTS == 2 else '2 (전체)' if N_DOTS == 6 else '?'}] "
      f"N_dots={N_DOTS}, channel_len={channel_len}nm, domain_l={domain_l}nm")

# Masks ---------------------------------------------------------------------------
channel = Polygon.box(channel_w, domain_l, name="channel").centered()
channel_mask = Mask("channel")
channel_mask.add_shape(channel)

oxide = Polygon.box(domain_w, domain_l, name="oxide").centered()
oxide_mask = Mask("oxide")
oxide_mask.add_shape(oxide)

source = (
    Polygon.box(channel_w, source_drain_w, name="source")
    .centered()
    .translated(0, -(domain_l - source_drain_w) / 2)
)
drain = (
    Polygon.box(channel_w, source_drain_w, name="drain")
    .centered()
    .translated(0, (domain_l - source_drain_w) / 2)
)
sd_mask = Mask("source_drain")
sd_mask.add_shapes([source, drain])

# --- Gate 좌표 계산 (B0 P1 B1 P2 ... P_N B_N 순서로 y_cursor를 전진시키며 배치) ---
y_cursor = -channel_len / 2
barrier_centers = []
plunger_centers = []
for i in range(N_DOTS):
    barrier_centers.append(y_cursor + barrier_w / 2)
    y_cursor += barrier_w + pitch
    plunger_centers.append(y_cursor + plunger_w / 2)
    y_cursor += plunger_w + pitch
barrier_centers.append(y_cursor + barrier_w / 2)  # 마지막 barrier B_{N_dots}
y_cursor += barrier_w
assert abs(y_cursor - channel_len / 2) < 1e-9, f"y_cursor={y_cursor}, expected {channel_len/2}"

gates = []
dots = []
for i, yc in enumerate(barrier_centers):
    gates.append(Polygon.box(domain_w, barrier_w, name=f"B{i}").centered().translated(0, yc))
for i, yc in enumerate(plunger_centers):
    gates.append(Polygon.box(domain_w, plunger_w, name=f"P{i+1}").centered().translated(0, yc))
    dots.append(Polygon.box(channel_w + pitch, QD_w, name=f"QD{i+1}").centered().translated(0, yc))

gate_mask = Mask("gates")
gate_mask.add_shapes(gates)
dot_mask = Mask("dots")
dot_mask.add_shapes(dots)

# Setup builder ---------------------------------------------------------------
# [WARNING] M12/M13/M17에서 확인된 바와 같이 Builder.view()/.view_shapes()는 이
# 환경에서 15분+/1100+ CPU초 걸려(메모리 문제 아님, VTK 오프스크린 렌더 병목)
# 사용하지 않는다. get_groups(sync=True)로 동기화만 수행.
builder = Builder()
builder.set_mesh_size(char_len).number_hull_surfaces(True)

builder.add_mask(channel_mask).add_mask(oxide_mask).add_mask(sd_mask).add_mask(
    gate_mask
).add_mask(dot_mask)
builder.use_all_masks()

# Build the BOX ---------------------------------------------------------------
builder.use_mask("oxide")
builder.set_z(-box_thick - channel_thick).extrude(box_thick)
builder.get_groups(sync=True)
builder.rename_group("oxide_bottom", "back_gate_bnd")
builder.dissolve_physical_group(lambda g: "oxide" in g.name and g.dim == 2)
builder.get_groups(sync=True)

# Build the channel -----------------------------------------------------------
builder.use_mask("channel")
builder.set_z(-channel_thick).extrude(channel_thick)
builder.dissolve_physical_group(lambda g: "channel" in g.name and g.dim == 2)
builder.get_groups(sync=True)

# Build the STI regions -------------------------------------------------------
builder.fill_mode()
builder.use_mask("oxide").set_z_from_group("channel", bottom=True)
builder.extrude(channel_thick)
builder.get_groups(sync=True)
builder.dissolve_physical_group(lambda g: "oxide" in g.name and g.dim == 2)
builder.dissolve_physical_group(lambda g: "channel" in g.name and g.dim == 2)
builder.get_groups(sync=True)

# Build source/drain (SD1/SD2 placeholder) regions -----------------------------
builder.displace_mode()
builder.set_z_from_group("channel", bottom=True)
builder.use_mask("source_drain")
builder.extrude(channel_thick)
builder.get_groups(sync=True)
builder.rename_group("source_side[3]", "source_bnd")
builder.rename_group("drain_side[1]", "drain_bnd")
builder.dissolve_physical_group(
    lambda g: "source" in g.name and g.dim == 2 and "bnd" not in g.name
)
builder.dissolve_physical_group(
    lambda g: "drain" in g.name and g.dim == 2 and "bnd" not in g.name
)
builder.get_groups(sync=True)
builder.number_hull_surfaces(False)

# Build gate-oxide layer ------------------------------------------------------
builder.use_mask("oxide").extrude(EOT_thick)
builder.dissolve_physical_group(lambda g: "oxide" in g.name and g.dim == 2)

# Deposit gates ---------------------------------------------------------------
builder.use_mask("gates")
builder.add_surface()
builder.get_groups(sync=True)

# Build dot regions -----------------------------------------------------------
builder.overlay_mode()
builder.use_mask("dots").set_z_from_group("channel", bottom=True)
builder.extrude(channel_thick + EOT_thick)
builder.get_groups(sync=True)

# Merge gate groups (B0..B_{N_dots}, P1..P_{N_dots}) ---------------------------
# [BASELINE] M17 1-fdsoi_builder.py와 동일한 substring-predicate 패턴을 그대로
# 재사용 (검증된 패턴, N_dots<=9에서 "B1"/"P1" 등이 다른 이름과 우연히 겹칠
# 위험 없음).
gate_names = [f"B{i}" for i in range(N_DOTS + 1)] + [f"P{i+1}" for i in range(N_DOTS)]
for gname in gate_names:
    target = f"barrier_gate_{gname[1:]}_bnd" if gname.startswith("B") else f"plunger_gate_{gname[1:]}_bnd"
    builder.merge_groups(lambda g, gname=gname: gname in g.name and g.dim == 2, target)

# [WARNING] M17에서 확인된 버그 패턴: naive "QD" predicate가 merge_groups가 이미
# 소비한 복합 이름(예: "B3.QD4_top")과 겹칠 수 있다 — 게이트 이름이 포함된 것은
# 명시적으로 제외해 순수 "QD_i_top/bottom" 잔여 surface만 dissolve한다.
builder.dissolve_physical_group(
    lambda g: "QD" in g.name
    and g.dim == 2
    and not any(gn in g.name for gn in gate_names)
)
builder.get_groups(sync=True)

# Save the mesh -----------------------------------------------------------------
# [WARNING] 한글 경로에서 gmsh XAO writer 실패 (M12~M17과 동일 버그) — C:\temp
# ASCII 스테이징 경로를 거쳐 프로젝트 폴더로 복사.
builder.mesh()
mesh_dir_ascii = Path(r"C:\temp\sixqubit")
mesh_dir_ascii.mkdir(parents=True, exist_ok=True)
tag = f"n{N_DOTS}"
builder.write(mesh_dir_ascii / f"array_{tag}.msh").write(mesh_dir_ascii / f"array_{tag}.xao")

out_meshes = script_dir / "meshes"
out_meshes.mkdir(parents=True, exist_ok=True)
shutil.copy(mesh_dir_ascii / f"array_{tag}.msh", out_meshes / f"array_{tag}.msh")
shutil.copy(mesh_dir_ascii / f"array_{tag}.xao", out_meshes / f"array_{tag}.xao")
print(f"Saved: {out_meshes / f'array_{tag}.msh'}")
