"""M07 Study 공통 모듈: Si/SiGe 원자 구조 생성과 Ge 판별 도우미.

02_tb_atoms_distribution.py, 03_tb_seed_ensemble.py가 함께 쓴다.

계면을 만드는 두 가지 방식
------------------------------
- "rough"  : 튜토리얼(SiGe Part 2) 방식. 우물 상·하부 계면을 RoughSurface로
             거칠게 만들고, 그 바깥 층을 일정 조성(Si 1-x, Ge x)의 무작위 합금으로 채운다.
             계면 자체는 원자 한 층 안에서 Si→SiGe로 급격히 바뀐다.
- "graded" : Atoms 5(spherical_dot) 방식을 층 구조에 적용. Ge 농도를 z의 함수
             x_Ge(z)로 주고, 원자 하나하나를 그 확률로 Si 또는 Ge로 뽑는다.
             계면 폭(상호확산) 효과를 볼 수 있다. 계면 거칠기는 넣지 않는다.

계면 프로파일 정의 (01 MVEMT 스크립트와 동일)
------------------------------
s(u) = 0.5 * (1 + tanh(u / w)),  w = INTERFACE_WIDTH (특성 길이)
우물 창(window) = s(z_top - z) * s(z - z_bot)      (우물 안 1, 밖 0)
x_Ge(z) = x_Ge_barrier * (1 - window(z))
→ 10 %→90 % 전이 폭 = 2·atanh(0.8)·w ≈ 2.2 w.  w = 0 이면 계단 함수.
"""

import numpy as np
from qtcad.atoms import Atoms
from qtcad.atoms.rough_surface import RoughSurface


# ----------------------------------------------------------------------------
# 계면 프로파일
# ----------------------------------------------------------------------------
def sigmoid(u, w):
    """0→1 전이 함수. w <= 0이면 계단 함수(u=0에서 0.5)."""
    if w <= 0:
        return np.heaviside(u, 0.5)
    return 0.5 * (1.0 + np.tanh(u / w))


def well_window(z, z_top, z_bot, w):
    """우물 안에서 1, 바깥에서 0인 창 함수 (단위는 입력과 동일)."""
    return sigmoid(z_top - z, w) * sigmoid(z - z_bot, w)


def ge_fraction_profile(z, z_top, z_bot, w, x_ge):
    """목표 Ge 조성 x_Ge(z)."""
    return x_ge * (1.0 - well_window(z, z_top, z_bot, w))


# ----------------------------------------------------------------------------
# 원자 구조 생성
# ----------------------------------------------------------------------------
def build_atoms(x, y, z, well_top, well_bot, x_ge, mode="rough",
                interface_width=0.0, hurst=0.3, rms=5e-10,
                seed=None, rough_seeds=(0, 1)):
    """Si/SiGe 우물 원자 구조를 만든다. 모든 길이는 m 단위.

    Args:
        x, y, z (np.ndarray): 원자 구조 상자의 [min, max] (m).
        well_top, well_bot (float): Si 우물 상·하부 계면의 평균 높이 (m).
        x_ge (float): 장벽(SiGe)의 Ge 조성 (예: 0.3).
        mode (str): "rough" 또는 "graded".
        interface_width (float): graded 모드의 계면 특성 길이 w (m).
        hurst, rms: rough 모드의 계면 거칠기 파라미터.
        seed (int | None): 합금 배치 난수 시드 (None이면 OS 엔트로피).
        rough_seeds (tuple): 상·하부 거친 계면의 난수 시드.

    Returns:
        atoms (Atoms), surfaces (dict): rough 모드면 {"top": ..., "bot": ...}.
    """
    atoms = Atoms(x=x, y=y, z=z, rng_seed=seed)
    surfaces = {}

    if mode == "rough":
        # 튜토리얼과 같은 순서: 거친 계면 생성 → 층 채우기
        surfaces["top"] = RoughSurface(hurst=hurst, mean=well_top, rms=rms,
                                       rng_seed=rough_seeds[0])
        surfaces["bot"] = RoughSurface(hurst=hurst, mean=well_bot, rms=rms,
                                       rng_seed=rough_seeds[1])
        alloy = np.array([["Si", 1.0 - x_ge], ["Ge", x_ge]])
        atoms.new_layer(z_top=z[1], z_bot=surfaces["top"], atom_species=alloy)  # spacer
        atoms.new_layer(z_top=surfaces["bot"], z_bot=z[0], atom_species=alloy)  # buffer
        # 두 층 사이(우물)는 기본값인 순수 Si로 남는다.

    elif mode == "graded":
        def ge_conc(xx, yy, zz):
            return float(ge_fraction_profile(zz, well_top, well_bot,
                                             interface_width, x_ge))

        def si_conc(xx, yy, zz):
            return 1.0 - ge_conc(xx, yy, zz)

        # 행마다 [원자 종, 농도 함수]. 위치를 주지 않으면 구조 전체에 적용된다.
        species = np.array([["Si", si_conc], ["Ge", ge_conc]], dtype=object)
        atoms.new_region(atom_species=species)

    else:
        raise ValueError(f"mode는 'rough' 또는 'graded'여야 합니다: {mode}")

    return atoms, surfaces


# ----------------------------------------------------------------------------
# 분석 도우미
# ----------------------------------------------------------------------------
def is_ge(species):
    """atom_species 배열에서 Ge인 원자를 True로 표시.
    문자열("Ge")과 원자번호(32) 표기를 모두 처리한다."""
    sp = np.asarray(species)
    if sp.dtype.kind in "iu":
        return sp == 32
    return np.char.strip(sp.astype(str)) == "Ge"


def positions_nm(atoms):
    """atom_pos를 nm 단위로 반환 (QTCAD Atoms의 좌표는 m 단위로 가정, 확인 출력 포함)."""
    pos = np.asarray(atoms.atom_pos, dtype=float)
    span = np.abs(pos).max()
    if span > 1e-6:  # m 단위라면 수십 nm = 수e-8 수준이어야 함
        raise RuntimeError(
            f"atom_pos의 최대 절댓값이 {span:.3e} 입니다. m 단위가 아닌 것 같으니 "
            "단위를 확인하고 positions_nm()의 변환을 고치세요.")
    return pos * 1e9
