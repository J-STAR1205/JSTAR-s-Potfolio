// Author: Kyrll Chudomirovic Flins
// SiO2-SiGe-Si-SiGe heterostructure
// Remark: length units are nanometers

SetFactory("Built-in");

Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// in-plane refinement
alpha  = 10;
fine   = 10 / alpha;   // 1.0 nm : 양자점 영역
medium = 20 / alpha;   // 2.0 nm : 게이트 주변
coarse = 40 / alpha;   // 4.0 nm : 외곽

// vertical refinement (분리)
beta = 1;
layers_ox         = beta * 9;    // 1.67 nm/층
layers_cap        = beta * 15;   // 1.87 nm/층
layers_well_above = beta * 5;    // 0.40 nm/층
layers_well       = beta * 20;   // 0.40 nm/층
layers_well_below = beta * 5;    // 0.40 nm/층
layers_buffer     = beta * 15;   // 6.7 nm/층

// lengthscales
size = 200;              // 유지
d = 40;                  // 논문 미기재 → 팀 마스크 설계의 플런저 폭 사용
gap = 5;                 // 논문 미기재 → 팀 마스크 설계의 게이트 간격 사용
thick_ox = 15;           // 아래 설명 참고 (5 또는 10)
thick_well_above = 2;    // 유지 (Ge 30% 동일)
thick_cap = 30 - thick_well_above;      // 28 nm, 논문 스페이서 30 nm와 일치
thick_well = 8;          // 논문: 8 nm 28Si
thick_well_below = 2;    // 유지
thick_buffer = 100 - thick_well_below;  // 98 nm

// plunger gate
Point(1) = {0, 0, 0, fine};
Point(2) = {20, 0, 0, fine};   // +x 방향, 반지름만큼
Point(3) = {0, 50, 0, fine};   // +y 방향
Point(4) = {-20, 0, 0, fine};   // -x 방향
Point(5) = {0, -50, 0, fine};   // -y 방향
Ellipse(1) = {2, 1, 3, 3};   // 시작:점2, 중심:점1, 장축기준점:점3, 끝:점3
Ellipse(2) = {3, 1, 3, 4};   // 시작:점3, 중심:점1, 장축기준점:점3, 끝:점4
Ellipse(3) = {4, 1, 3, 5};   // 시작:점4, 중심:점1, 장축기준점:점3, 끝:점5
Ellipse(4) = {5, 1, 3, 2};   // 시작:점5, 중심:점1, 장축기준점:점3, 끝:점2
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// annular gap between plunger gate and confinement gate
Point(6) = {25, 0, 0, medium};    // +x, 바깥 타원
Point(7) = {0, 55, 0, medium};    // +y, 바깥 타원
Point(8) = {-25, 0, 0, medium};    // -x
Point(9) = {0, -55, 0, medium};    // -y
Ellipse(5) = {6, 1, 7, 7};         // D단계와 같은 패턴 (장축기준점=점7)
Ellipse(6) = {7, 1, 7, 8};
Ellipse(7) = {8, 1, 7, 9};
Ellipse(8) = {9, 1, 7, 6};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2, 1};         // 바깥루프(2)에서 안쪽루프(1)를 구멍으로 뺌

// confinement gate
Point(10) = {-size/2, size/2, 0, coarse};   // 좌상
Point(11) = {size/2, size/2, 0, coarse};    // 우상
Point(12) = {size/2, -size/2, 0, coarse};   // 우하
Point(13) = {-size/2, -size/2, 0, coarse};  // 좌하
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};
Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3, 2};   // 바깥루프(3)에서 안쪽루프(2, 갭 바깥경계)를 구멍으로 뺌

// extrusions along z axis
thick_list = {thick_ox, thick_cap, thick_well_above, thick_well, thick_well_below, thick_buffer};    // 힌트: B/C단계 변수 이름을 순서대로 (산화막→캡→우물위→우물→우물아래→버퍼)
layers_list = {layers_ox, layers_cap, layers_well_above, layers_well, layers_well_below, layers_buffer};   // 힌트: 위와 같은 순서로 레이어수 변수들

// first extrusions
extrusion_01[] = Extrude {0, 0, -thick_list[0]} {Surface{1}; Layers{layers_list[0]};};   // D의 면 태그
extrusion_02[] = Extrude {0, 0, -thick_list[0]} {Surface{2}; Layers{layers_list[0]};};   // E의 면 태그
extrusion_03[] = Extrude {0, 0, -thick_list[0]} {Surface{3}; Layers{layers_list[0]};};   // F의 면 태그

// remaining extrusions
For id In {1:5}
    extrusion_01[] = Extrude {0, 0, -thick_list[id]} {Surface{extrusion_01[0]}; Layers{layers_list[id]};};
    extrusion_02[] = Extrude {0, 0, -thick_list[id]} {Surface{extrusion_02[0]}; Layers{layers_list[id]};};
    extrusion_03[] = Extrude {0, 0, -thick_list[id]} {Surface{extrusion_03[0]}; Layers{layers_list[id]};};
EndFor

Physical Surface("gate_plunger") = {1};       // D단계 면 태그
Physical Surface("gate_confinement") = {3};   // F단계 면 태그

Physical Volume("oxide") = {1, 2, 3};          // 층1 전체 3기둥
Physical Volume("cap_top") = {4, 5, 6};        // 층2(cap) 전체 3기둥
Physical Volume("cap_qd") = {7};               // 층3(well_above) D기둥 ← QD 서브메시용
Physical Volume("cap_bot") = {8, 9};           // 층3의 E+F기둥
Physical Volume("well_qd") = {10};             // 층4(진짜 우물) D기둥 ← 핵심 QD
Physical Volume("well") = {11, 12};            // 층4의 E+F기둥
Physical Volume("buffer_qd") = {13};           // 층5(well_below) D기둥 ← QD 서브메시용
Physical Volume("buffer_top") = {14, 15};      // 층5의 E+F기둥
Physical Volume("buffer_bot") = {16, 17, 18};  // 층6(buffer) 전체 3기둥

Mesh 1;
Mesh 2;
Mesh 3;
Save "multiscale.msh";
