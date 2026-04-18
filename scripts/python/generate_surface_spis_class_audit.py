from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SPIS_ROOT = REPO_ROOT / "ref" / "SPIS" / "num-master" / "src" / "main" / "java" / "spis"
DOCS_DIR = REPO_ROOT / "Toolkit" / "Surface Charging" / "docs"
CLASS_AUDIT_JSON = DOCS_DIR / "surface_spis_class_audit.json"
CLASS_AUDIT_CSV = DOCS_DIR / "surface_spis_class_audit.csv"
CLASS_AUDIT_SUMMARY = DOCS_DIR / "surface_spis_class_audit_summary.md"

FAMILY_MATRIX_PATH = REPO_ROOT / "scripts" / "run" / "surface_spis_family_coverage_matrix.json"
SEALOFF_PATH = REPO_ROOT / "build_codex" / "surface_round47_sealoff" / "surface_round26_sealoff.json"
PLAN_DOC_PATH = REPO_ROOT / "Toolkit" / "Surface Charging" / "docs" / "surface_num_alignment_refactor_plan.md"
CONTRACT_PATH = REPO_ROOT / "documentation" / "contracts" / "surface_family_coverage_matrix_v1.json"
GATE_PATH = REPO_ROOT / "build_codex" / "surface_family_coverage_gate_round47" / "surface_family_coverage_gate.json"
REFERENCE_MATRIX_PATH = REPO_ROOT / "scripts" / "run" / "surface_reference_matrix.json"

STATUS_EQUIVALENT = "已等价实现"
STATUS_AGGREGATED = "聚合等价实现"
STATUS_APPROX = "近似实现"
STATUS_MISSING = "尚无实现"
STATUS_OUT = "不纳入当前目标范围"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def norm_path_string(path: str) -> str:
    return path.replace("\\", "/")


FAMILY_RULES = [
    ("Surf/Material/", "surf_material", "Surf/Material"),
    ("Surf/SurfInteract/", "surf_interact", "Surf/SurfInteract"),
    ("Surf/SurfDistrib/", "surf_distrib", "Surf/SurfDistrib"),
    ("Surf/SurfField/", "surf_field", "Surf/SurfField"),
    ("Surf/SurfMesh/", "surf_mesh", "Surf/SurfMesh"),
    ("Circ/CircField/", "circ_circfield", "Circ/CircField"),
    ("Circ/DIDV/", "circ_didv", "Circ/DIDV"),
    ("Circ/Circ/", "circ_circuit", "Circ/Circ"),
    ("Solver/Circuit/", "solver_circuit", "Solver/Circuit"),
    ("Solver/ElectroMag/", "solver_electromag", "Solver/ElectroMag"),
    ("Solver/Matter/", "solver_matter", "Solver/Matter"),
    ("Solver/Util/", "solver_util", "Solver/Util"),
    ("Top/Transition/", "top_transition", "Top/Transition"),
    ("Top/Simulation/", "top_simulation", "Top/Simulation"),
    ("Top/Plasma/", "top_plasma", "Top/Plasma"),
    ("Top/SC/", "top_sc", "Top/SC"),
    ("Top/Grid/", "top_grid", "Top/Grid"),
    ("Top/Default/", "top_default", "Top/Default"),
    ("Top/Top/", "top_top", "Top/Top"),
    ("Util/DistribFunc/", "util_distrib_func", "Util/DistribFunc"),
    ("Util/Exception/", "util_exception", "Util/Exception"),
    ("Util/Instrument/", "util_instrument", "Util/Instrument"),
    ("Util/Func/", "util_func", "Util/Func"),
    ("Util/Part/", "util_part", "Util/Part"),
    ("Util/Phys/", "util_phys", "Util/Phys"),
    ("Util/Sampler/", "util_sampler", "Util/Sampler"),
    ("Util/Monitor/", "util_monitor", "Util/Monitor"),
    ("Util/io/", "util_io", "Util/io"),
    ("Util/Table/", "util_table", "Util/Table"),
    ("Util/Matrix/", "util_matrix", "Util/Matrix"),
    ("Util/Vect/", "util_vect", "Util/Vect"),
    ("Util/List/", "util_list", "Util/List"),
    ("Util/OcTreeMesh/", "util_octree_mesh", "Util/OcTreeMesh"),
    ("Util/OcTree/", "util_octree", "Util/OcTree"),
    ("Vol/BC/", "vol_bc_abstract", "Vol/BC(Abstract)"),
    ("Vol/VolField/", "vol_field_abstract", "Vol/VolField(Abstract)"),
    ("Vol/VolDistrib/", "vol_distrib_long_tail", "Vol/VolDistrib(LongTail)"),
    ("Vol/VolInteract/", "vol_interact_long_tail", "Vol/VolInteract(LongTail)"),
    ("Vol/Geom/", "vol_geom", "Vol/Geom"),
    ("Vol/VolMesh/", "vol_mesh", "Vol/VolMesh"),
    ("Vol/", "vol_full_stack", "Vol/BC+Field+Distrib+Interact"),
]


FAMILY_DEFAULTS = {
    "osgi_out_of_scope": {
        "default_status": STATUS_OUT,
        "mapping": "OSGi bundle 壳层；不属于当前 Surface Charging 主线数值对齐对象。",
        "difference": "该类属于 SPIS 框架/打包壳层，不在当前数值主线审计范围。",
        "next_action": "保持排除在当前审计范围之外。",
    },
    "surf_material": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Tools/Material 中的材料模型与参数资产 + DensePlasmaSurfaceCharging 中的材料路由。",
        "difference": "当前项目按 C++ 材料模型与参数资产组织，不保留 SPIS Java 单类层级。",
        "next_action": "保持材料参数资产、单测和 reference matrix 常绿。",
    },
    "surf_interact": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Tools/Material 的 SurfaceInteraction 家族 + DensePlasmaSurfaceCharging 的 interactor 路由。",
        "difference": "当前项目把交互器与 yield/helper 函数拆到 C++ 内核与 glue 逻辑，不保持 Java 插件对象树。",
        "next_action": "继续维护显式 family route 与 Material/Smoke 证据链。",
    },
    "surf_distrib": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Tools/Particle 的 SurfaceDistributionFunction 家族 + DensePlasmaSurfaceCharging 分布路由。",
        "difference": "当前项目将 SurfDistrib 多模型收口到统一分布函数实现与运行时选择器中。",
        "next_action": "保留 distribution family 的 smoke、reference matrix 与参数兼容性。",
    },
    "surf_field": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Tools/FieldSolver 的 SurfaceBarrierModels/BarrierCurrentScaler + DensePlasmaSurfaceCharging 元数据导出。",
        "difference": "当前项目以 scaler 组合内核承载 SPIS SurfField 家族，不保留 Java 抽象层级外形。",
        "next_action": "继续维护 scaler family 的 direct route 与 barrier smoke。",
    },
    "surf_mesh": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceChargingKernelFramework + DensePlasmaSurfaceCharging + Main 配置入口中的 surface mesh/runtime 组织层。",
        "difference": "当前项目把 SurfMesh 语义融入场景装配与 mesh/runtime glue，而非保留独立 Java mesh 对象树。",
        "next_action": "保持 structured-topology route、sidecar 与 graph/report smoke。",
    },
    "circ_circuit": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Tools/Coupling 的 SurfaceCircuitCoupling/CircuitAssembly + DensePlasmaSurfaceCharging 的电路装配。",
        "difference": "当前项目以通用耦合内核承载电路元件和激励，不保持 SPIS 电路类的逐对象封装外形。",
        "next_action": "继续维护 multi-patch circuit 场景与 family metadata。",
    },
    "circ_circfield": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceCircuitCoupling + SurfaceCurrentLinearizationKernel 中的 CircField/CurrentDistribution 视图。",
        "difference": "当前项目把 CircField 语义并入线性化与 current distribution 内核，而非独立 Java field 类树。",
        "next_action": "保持 circuit field route 与 finite-difference 回归。",
    },
    "circ_didv": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceCircuitCoupling 中的 DIDV 组合、约化与矩阵构造路径。",
        "difference": "当前项目用组合器和 reduction kernel 实现 DIDV，不保留所有 Java DIDV 包装类的形状。",
        "next_action": "继续维护 DIDV reduction 与 reference matrix 对照。",
    },
    "solver_circuit": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceSolverFacade 中的 circuit solve 视图 + DensePlasmaSurfaceCharging runtime solver route。",
        "difference": "当前项目以 facade + runtime route 承载 circuit solver 语义，而不是单独的 Java solver 类。",
        "next_action": "保持 solver family metadata 与 circuit 路径 smoke。",
    },
    "solver_util": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Tools/Solver LinearSystem + Tools/FieldSolver PoissonSolver + PICCycle 中的共用求解工具路径。",
        "difference": "当前项目把 Solver/Util 能力拆到线性系统、泊松求解和 PIC 公共流程，不保留 Java util 类型层级。",
        "next_action": "保持 least-square/linear system 相关单测与 solver metadata。",
    },
    "solver_electromag": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceSolverFacade + SurfaceChargingKernelFramework 中的电磁/泊松求解路由。",
        "difference": "当前项目通过 facade 和线性求解策略承载电磁求解，不维持 Java solver 继承树镜像。",
        "next_action": "继续维护 Poisson/electromag route 与 reference matrix 阈值。",
    },
    "solver_matter": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceSolverFacade + DensePlasmaSurfaceCharging 中的 particle push/runtime matter route。",
        "difference": "当前项目以 runtime matter route 承载粒子推进家族，不保留 Java pusher 子类的逐类外形。",
        "next_action": "继续维护 live-PIC / hybrid matter 路由与 smoke 证据。",
    },
    "top_transition": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceTransitionEngine 中的 transition/updater/lifecycle 对象层 + DensePlasmaSurfaceCharging 元数据导出。",
        "difference": "当前项目已闭环 transition 语义，但对象层级被收口为统一 transition engine。",
        "next_action": "保持 transition lifecycle smoke 与 object-layer 测试。",
    },
    "top_top": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceScenarioCatalog、SurfaceSimulationRunner、SurfaceRuntimePlan、Main 入口中的 top-level 场景组织层。",
        "difference": "当前项目以场景目录和 runtime plan 替代 SPIS 顶层 UI/top 对象树。",
        "next_action": "继续维护 top route evidence 与 scenario/runtime 入口契约。",
    },
    "top_simulation": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceSimulationRunner + SurfaceRuntimePlan + DensePlasmaSurfaceCharging 模拟主线。",
        "difference": "当前项目将 Simulation 家族压缩为统一 simulation runner/runtime plan，不保留 Java listener/init 对象层。",
        "next_action": "保持 simulation route、reference matrix 与 lifecycle hooks 回归。",
    },
    "top_plasma": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "DensePlasmaSurfaceCharging + SurfaceRuntimePlan + Plasma Analysis 配置对象中的 plasma/environment 路由。",
        "difference": "当前项目按环境参数与 runtime 计划承载 plasma 家族，而不是逐 Java plasma 类镜像。",
        "next_action": "继续维护 GEO/LEO plasma case 与 distribution/runtime 兼容性。",
    },
    "top_sc": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceCircuitCoupling + DensePlasmaSurfaceCharging 中的 spacecraft/body/patch/interface 组织层。",
        "difference": "当前项目把 SC 家族融入 circuit/body/patch 运行时装配，不保留 Java SC 类树外形。",
        "next_action": "保持 multi-body / multi-patch 场景和 SC metadata。",
    },
    "top_grid": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Main/SurfaceScenarioLoader + SurfaceFieldVolumeBridge + DensePlasmaSurfaceCharging 中的 grid/boundary mapping。",
        "difference": "当前项目把 grid 组织层并入 loader 与 bridge，不保留 Java grid 对象层。",
        "next_action": "继续维护 boundary mapping、grid metadata 和 mesh/volume 路由。",
    },
    "top_default": {
        "default_status": STATUS_APPROX,
        "mapping": "SurfaceScenarioCatalog + SurfaceRuntimePlan + DensePlasmaSurfaceCharging 中的通用对象、参数与默认值组织层。",
        "difference": "Top/Default 的大量框架型 Java 基类被当前项目折叠为配置字段、runtime structs 与 glue 逻辑，类级镜像最弱。",
        "next_action": "若要提升 class 级可信度，优先补 Top/Default 对象层索引与显式职责映射。",
    },
    "util_full_stack": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Tools/Basic、Geometry、Diagnostics、Output、Particle、Solver、Mesh 中的基础设施集合。",
        "difference": "当前项目已按公共库方式承载 Util/**，但不保留 SPIS Java util 包的逐类对象组织。",
        "next_action": "将 class-level 审计结果作为 util 子包持续拆分索引基线。",
    },
    "util_distrib_func": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceDistributionFunction 与相关分布构造工具。",
        "difference": "当前项目以统一分布函数内核承载多个 Java DF 类，逐类对象镜像较弱。",
        "next_action": "维持 distribution family 测试与 pressure/monitor 相关语义说明。",
    },
    "util_exception": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "ErrorCode、ErrorHandling、VoidResult 及 DensePlasmaSurfaceCharging 错误桥接。",
        "difference": "异常语义已覆盖，但异常类型被收口为 C++ 错误码/异常辅助设施。",
        "next_action": "保持异常路径测试与错误码契约稳定。",
    },
    "util_instrument": {
        "default_status": STATUS_APPROX,
        "mapping": "DensePlasmaSurfaceCharging 的 observer/instrument metadata、Main 配置入口与 smoke 断言。",
        "difference": "当前项目已具备 instrument family 运行时证据，但多数 SPIS instrument 单类未下沉为独立 C++ 类型。",
        "next_action": "如需提升类级对齐，优先拆 observer/instrument catalog 与专门单测。",
    },
    "util_func": {
        "default_status": STATUS_APPROX,
        "mapping": "MathUtils、NumericAggregation、StringTokenUtils、ModelRegistry 等基础函数与工具库。",
        "difference": "Util/Func 在当前项目中主要被折叠为基础函数库和公式 helpers，未保留逐 Java 函数对象类层级。",
        "next_action": "如需提高 class 级可信度，优先为物理相关函数建立显式索引与命名桥接。",
    },
    "util_part": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "ParticleDefinitions、ParticleSource 与 DensePlasmaSurfaceCharging 中的粒子数据组织层。",
        "difference": "当前项目已承载粒子数据与 source family，但 Java util part 类未逐个镜像为独立 C++ 类型。",
        "next_action": "保持粒子定义单测与 source/runtime 兼容性。",
    },
    "util_phys": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Constants、MathUtils、ParticleDefinitions 中的物理常数与物理 helper。",
        "difference": "当前项目把 Phys 家族收口到基础常数与数学/粒子 helper 中，而非独立 Java phys 类。",
        "next_action": "维持物理常数与基础公式的稳定接口。",
    },
    "util_sampler": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceDistributionFunction + PicMccSurfaceCurrentSampler + DensePlasmaSurfaceCharging 采样路径。",
        "difference": "当前项目具备 sampler 语义闭环，但采样器类被整合为统一 surface/PIC 采样内核。",
        "next_action": "继续维护 sampler metadata 与 sampling smoke。",
    },
    "util_monitor": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Monitor、FieldMonitor、MonitorManager 与 DensePlasmaSurfaceCharging metadata。",
        "difference": "当前项目已实现监测能力，但 monitor class hierarchy 以通用监测器管理器承载。",
        "next_action": "保持 Diagnostics 单测与 monitor sidecar 契约。",
    },
    "util_io": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "ResultExporter、CSVExporter、SurfaceScenarioLoader、DensePlasmaSurfaceCharging 输出与输入桥接。",
        "difference": "当前项目以统一 IO/exporter/loader 流程承载 Util/io，不保留所有 Java 读写器类。",
        "next_action": "保持 IO 输出契约、CLI loader 与 reference matrix 兼容性。",
    },
    "util_table": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "ParticleSource、SurfaceDistributionFunction、DataAnalyzer 与 DensePlasmaSurfaceCharging 的表格化数据路径。",
        "difference": "当前项目把表格/插值能力并入 source/distribution/output 组件，不保留 Java table 对象层。",
        "next_action": "继续维护表格类工具的单测与数据输出一致性。",
    },
    "util_matrix": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "LinearSystem、SparseMatrix、BandMatrix、SurfaceSolverFacade 中的矩阵与线性代数内核。",
        "difference": "当前项目已具备矩阵能力，但 Java matrix class 家族被折叠为求解器公共数据结构。",
        "next_action": "保持矩阵/线性系统测试与 solver 路由稳定。",
    },
    "util_vect": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "Vector3D、Point3D 与 DensePlasmaSurfaceCharging 中的矢量计算路径。",
        "difference": "当前项目以 Geometry 基础类型统一承载矢量家族，而非逐 Java vect class 镜像。",
        "next_action": "保持 Geometry 单测与向量 API 稳定。",
    },
    "util_list": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "ParticleSource、SurfaceDistributionFunction 与 DensePlasmaSurfaceCharging 中的列表/集合式 source 组织层。",
        "difference": "当前项目将 List 家族折叠进粒子源与运行时容器，不保留 Java list util 类形态。",
        "next_action": "保持 ParticleSource 与 distribution 相关单测。",
    },
    "util_octree": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "MeshPartitioning 与 DensePlasmaSurfaceCharging 中的 octree 路由。",
        "difference": "当前项目以 mesh partitioning 内核承载 octree 语义，而不是保留 Java oc-tree 类型层级。",
        "next_action": "保持 mesh partitioning 测试和 octree metadata。",
    },
    "util_octree_mesh": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "MeshAlgorithms、MeshPartitioning 与 DensePlasmaSurfaceCharging 中的 octree-mesh 路由。",
        "difference": "当前项目以 mesh algorithm/partitioning 组合承载 octree mesh，不保留逐 Java mesh wrapper 类。",
        "next_action": "继续维护 Mesh 单测、octree mesh metadata 与相关 smoke。",
    },
    "vol_full_stack": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "SurfaceFieldVolumeContracts + SurfaceFieldVolumeBridge + DensePlasmaSurfaceCharging native-volume route。",
        "difference": "当前项目以 bridge/runtime parity 方式承载 Vol 全栈，而不是 SPIS 原生对象层直译。",
        "next_action": "把当前 class-level 审计作为后续 native-volume 对象层显式化基线。",
    },
    "vol_bc_abstract": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "BoundaryCondition、PoissonSolver、LinearSystem、SurfaceFieldVolumeBridge 中的 BC/VoltageGenerator 路径。",
        "difference": "当前项目已具备 BC family 语义与 override route，但仍以 bridge/runtime 对象承载，不保留 Java BC 包装层。",
        "next_action": "继续维护 native volume BC override、smoke 与 sealoff 基线。",
    },
    "vol_field_abstract": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "VolField、SurfaceFieldVolumeContracts、SurfaceFieldVolumeBridge 与 DensePlasmaSurfaceCharging 的 field route。",
        "difference": "当前项目已具备 VolField family 路由与 metadata，但类级对象层被合并到 field bridge 内核。",
        "next_action": "保持 native volume field override 与 bridge 测试。",
    },
    "vol_distrib_long_tail": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "VolDistrib、SurfaceFieldVolumeContracts、SurfaceFieldVolumeBridge、Main 配置入口与 DensePlasmaSurfaceCharging 分布路由。",
        "difference": "当前项目已闭环 VolDistrib 长尾 family，但多个 Java 分布类统一收口到配置+bridge route。",
        "next_action": "继续维护 native volume distribution override 与长尾 reference cases。",
    },
    "vol_interact_long_tail": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "VolInteract、SurfaceFieldVolumeBridge 与 DensePlasmaSurfaceCharging interaction route。",
        "difference": "当前项目以 interaction kernel + bridge 组合承载 VolInteract 长尾 family，不保留全部 Java interactor 外形。",
        "next_action": "保持 interaction long-tail override、FieldSolver 单测与 reference matrix。",
    },
    "vol_geom": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "MeshParsing、SurfaceFieldVolumeBridge 与 DensePlasmaSurfaceCharging 中的体域几何语义。",
        "difference": "当前项目把 Vol/Geom 语义并入 mesh parsing 与 bridge/runtime，而不是独立 Java geometry 对象层。",
        "next_action": "保持 native-volume geometry metadata 与 mesh route smoke。",
    },
    "vol_mesh": {
        "default_status": STATUS_AGGREGATED,
        "mapping": "MeshParsing、SurfaceFieldVolumeBridge 与 DensePlasmaSurfaceCharging 中的体域 mesh 路由。",
        "difference": "当前项目已具备 VolMesh 路由，但 mesh class 语义更多由 parsing/bridge 承载，而非 Java mesh 树镜像。",
        "next_action": "继续维护 volume mesh metadata、bridge 测试与 sealoff。",
    },
}


DIRECT_CLASS_OVERRIDES = {
    "Surf/Material/BasicMaterialModel.java": "BasicSurfaceMaterialModel + SurfaceMaterialModel family route。",
    "Surf/Material/ErosionMaterialModel.java": "ErosionSurfaceMaterialModel + erosion 参数资产与 material route。",
    "Surf/Material/ErosionParamSet.java": "ErosionParamSet 参数资产与材料参数载体。",
    "Surf/SurfMesh/ThreeDUnstructSurfMesh.java": "SurfaceChargingKernelFramework 中的 structured/unstructured surface mesh route。",
    "Circ/Circ/Circuit.java": "SurfaceCircuitCoupling 中的 CircuitAssembly/Circuit family 视图。",
    "Solver/Circuit/CircSolve.java": "SurfaceSolverFacade 中的 solver circuit family 视图。",
    "Solver/Util/LeastSquare.java": "LinearSystem/solver util 路径中的 least-square 支撑能力。",
    "Solver/Util/SolverUtil.java": "SurfaceSolverFacade + LinearSystem 的 solver util 路径。",
    "Solver/ElectroMag/Poisson/PoissonSolver.java": "SurfaceSolverFacade + PoissonSolver 求解路径。",
    "Solver/ElectroMag/Poisson/PotPoissonSolver.java": "SurfaceSolverFacade 中的 potential/Poisson route。",
    "Solver/ElectroMag/Poisson/ConjGrad3DUnstructPoissonSolver.java": "IterativeLinearSolvers/Poisson solver route。",
    "Solver/ElectroMag/PoissonMagnetostaticSolver.java": "SurfaceSolverFacade 中的 electromag family route。",
    "Vol/VolField/VolField.java": "VolField 抽象层与 native volume field route。",
    "Vol/VolInteract/VolInteractor.java": "VolInteract 抽象层与 interaction family route。",
}


DIRECT_BY_FAMILY = {
    "surf_material": {
        "BasicMaterialModel",
        "ErosionMaterialModel",
        "ErosionParamSet",
    },
    "surf_interact": {
        "MaterialInteractor",
        "MaterialDFInteractor",
        "GenericDFInteractor",
        "MaxwellianInteractor",
        "MaxwellianInteractorWithRecollection",
        "MultipleInteractor",
        "MultipleMaxwellianInteractor",
        "YieldInteractor",
        "ReflectionInteractor",
        "ErosionInteractor",
        "ImprovedPhotoEmInteractor",
        "TabulatedSEYModel",
        "RecollectedSEYModel",
        "DefaultPEEModel",
        "DefaultSEEEYModel",
        "DefaultSEEPModel",
        "DefaultErosionModel",
        "BasicInducedConductInteractor",
        "Device",
    },
    "surf_distrib": {
        "AxisymTabulatedVelocitySurfDistrib",
        "FluidSurfDistrib",
        "FowlerNordheimSurfDistrib",
        "GlobalMaxwellBoltzmannSurfDistrib",
        "GlobalMaxwellBoltzmannSurfDistrib2",
        "GlobalMaxwellSurfDistrib",
        "LocalMaxwellSurfDistrib",
        "LocalMaxwellSurfDistrib2",
        "LocalModifiedPearsonIVSurfDistrib",
        "LocalTabulatedSurfDistrib",
        "MaxwellianThruster",
        "MultipleSurfDistrib",
        "NonLocalizedSurfDistrib",
        "NonPICSurfDistrib",
        "PICSurfDistrib",
        "RecollMaxwellSurfDistrib",
        "TestableSurfDistrib",
        "TestableForASurfDistrib",
        "TwoAxesTabulatedVelocitySurfDistrib",
        "UniformVelocitySurfDistrib",
    },
    "surf_field": {
        "AutomaticBarrierCurrentScaler",
        "BarrierCurrentScaler",
        "CurrentScalerFromCurrentVariation",
        "CurrentScalerFromLocalCurrentVariation",
        "FowlerNordheimCurrentScaler",
        "GlobalTempCurrentScaler",
        "LTE_OML_CurrentScaler",
        "MultipleCurrentScaler",
        "OMLCurrentScaler",
        "SecondaryRecollectionCurrentScaler",
        "SmoothedGlobalTempCurrentScaler",
        "VariableBarrierCurrentScaler",
        "VariableBarrierCurrentScaler2",
    },
    "circ_circuit": {
        "Circuit",
        "RCCabsCirc",
        "RLCCirc",
        "SIN",
        "PULSE",
        "PWL",
        "EXP",
        "SurfaceComponent",
    },
    "circ_circfield": {
        "CircField",
        "CurrentDistribution",
        "DirCircField",
    },
    "circ_didv": {
        "DIDV",
        "DIDVOnCirc",
        "MultipleDIDV",
        "MultipleSurfDIDV",
        "ReducableDIDV",
        "SurfDIDV",
        "SurfDIDVFromCurrentVariation",
        "SurfDIDVFromLocalCurrentVariation",
        "SurfDIDVFromMatrices",
        "WeightedSurfDIDV",
    },
    "solver_circuit": {"CircSolve"},
    "solver_util": {"LeastSquare", "SolverUtil"},
    "solver_electromag": {
        "PoissonMagnetostaticSolver",
        "PoissonSolver",
        "PotPoissonSolver",
        "ConjGrad3DUnstructPoissonSolver",
    },
    "top_transition": {
        "BasicEclipseExit",
        "ConductivityEvolution",
        "Finalization",
        "LangmuirProbeTransition",
        "LocalTimeTransition",
        "RCCabsSCUpdater",
        "SheathOrPresheathPoissonBCUpdater",
        "SimulationParamUpdater",
        "SourceFluxUpdater",
        "SpinningSpacecraft",
        "SunFluxIntensityUpdater",
        "SunFluxUpdater",
        "TransientArtificialSources",
        "TransitionObserver",
    },
    "top_simulation": {
        "Simulation",
        "PlasmaScSimulation",
        "SimulationFromUIParams",
    },
    "vol_bc_abstract": {
        "DirichletPoissonBC",
        "MatterBC",
        "TestableMatterBC",
        "VoltageGenerator",
        "VoltageDependentMBC",
        "MixedDirichletFourierPoissonBC",
        "FourierPoissonBC",
        "SurfDistribMatterBC",
        "OneSurfDistribTestableMatterBC",
        "CapacitiveVoltageGenerator",
    },
    "vol_field_abstract": {
        "ScalVolField",
        "VectVolField",
        "EField",
        "DirVectVolField",
        "DirEField",
        "PotEField",
        "PotVectVolField",
        "BField",
        "UniformBField",
        "SolenoidBField",
        "DipolarBField",
        "MultipleEField",
        "EMField",
    },
    "vol_distrib_long_tail": {
        "PICVolDistrib",
        "PICVolDistribNoAcc",
        "PICVolDistribUpdatable",
        "SmartPICVolDistrib",
        "CompositeVolDistrib",
        "BackTrackingVolDistrib",
        "BacktrackingPICCompositeVolDistrib",
        "BacktrackingBoltzmannCompositeVolDistrib",
        "NonPICVolDistrib",
        "AnalyticVolDistrib",
        "MultipleVolDistrib",
        "LocalMaxwellVolDistrib",
        "LocalMaxwellBoltzmannVolDistrib",
        "GlobalMaxwellBoltzmannVolDistrib",
        "PICBoltzmannVolDistrib",
        "SteadyMaxwellBoltzmannVolDistrib",
        "UnlimitedGlobalMaxwellBoltzmannVolDistrib",
        "SurfaceLimitedGlobalMaxwellBoltzmannVolDistrib",
        "TrunckatedGlobalMaxwellBoltzmannVolDistrib",
        "ImplicitableVolDistrib",
        "Updatable",
    },
    "vol_interact_long_tail": {
        "VolInteractor",
        "VolInteractModel",
        "VolInteractionAlongTrajectory",
        "ElasticCollisions",
        "ConstantIonizationInteractor",
        "MCCInteractor",
        "CEXInteractor",
        "PhotoIonization",
        "TrajectoryInteractionFromField",
        "SpinningSpacecraftTrajectory",
    },
}


APPROX_BY_FAMILY = {
    "top_default": set(),
    "util_func": set(),
    "util_instrument": set(),
}


PATH_STATUS_OVERRIDES = {
    "osgi/SpisBundle.java": STATUS_OUT,
    "Top/Top/UIInvokable.java": STATUS_APPROX,
    "Top/Top/SpisTopMenu.java": STATUS_APPROX,
    "Top/Simulation/SimulationListener.java": STATUS_APPROX,
    "Top/Default/Common.java": STATUS_APPROX,
    "Top/Default/Credits.java": STATUS_APPROX,
    "Top/Default/Duplicable.java": STATUS_APPROX,
    "Top/Default/Listener.java": STATUS_APPROX,
    "Top/Default/ThreadManager.java": STATUS_APPROX,
    "Top/Default/SimulationObject.java": STATUS_APPROX,
    "Top/Default/NamedObject.java": STATUS_APPROX,
    "Top/Default/Parameter.java": STATUS_APPROX,
    "Top/Default/GlobalParameter.java": STATUS_APPROX,
    "Top/Default/LocalParameter.java": STATUS_APPROX,
    "Top/Default/InstrumentParameter.java": STATUS_APPROX,
    "Top/Default/PlasmaIO.java": STATUS_APPROX,
}


FAMILY_ATTENTION_NOTES = {
    "top_default": "Top/Default 在 family gate 上已封板，但 class 级对象层镜像最弱，建议单独补一份对象层职责索引。",
    "util_func": "Util/Func 大量 Java 函数对象已折叠为基础数学/物理 helper，若需要 class 级可信度，应补命名桥接与函数索引。",
    "util_instrument": "Util/Instrument 目前更接近 family route 完成，尚未形成逐 instrument 的 C++ 类型与专门单测矩阵。",
    "solver_matter": "Solver/Matter 已在 family 层闭环，但具体 pusher 子类仍以 runtime route 承载为主，类级显式对象层可继续增强。",
    "vol_full_stack": "Vol 全栈当前主要依赖 bridge/native-volume route 封板，若进一步追求 class 级 1:1，应继续拆 native volume 对象层。",
    "top_top": "Top/Top 的 UI/top-shell 类在当前项目中被 runtime plan 和 scenario runner 吸收，适合补一份 shell-to-runtime 对照。",
}


def resolve_family(spis_path: str) -> tuple[str, str]:
    if spis_path.startswith("osgi/"):
        return "osgi_out_of_scope", "osgi"
    for prefix, family_id, family_name in FAMILY_RULES:
        if spis_path.startswith(prefix):
            return family_id, family_name
    if spis_path.startswith("Util/"):
        return "util_full_stack", "Util/**"
    raise ValueError(f"Unable to resolve family for {spis_path}")


def infer_status(spis_path: str, family_id: str) -> str:
    class_name = Path(spis_path).stem
    if spis_path in PATH_STATUS_OVERRIDES:
        return PATH_STATUS_OVERRIDES[spis_path]
    if class_name in DIRECT_BY_FAMILY.get(family_id, set()):
        return STATUS_EQUIVALENT
    if family_id in APPROX_BY_FAMILY:
        return STATUS_APPROX
    return FAMILY_DEFAULTS.get(family_id, FAMILY_DEFAULTS["util_full_stack"])["default_status"]


def family_mapping_text(family_id: str, spis_path: str) -> str:
    if spis_path in DIRECT_CLASS_OVERRIDES:
        return DIRECT_CLASS_OVERRIDES[spis_path]
    return FAMILY_DEFAULTS.get(family_id, FAMILY_DEFAULTS["util_full_stack"])["mapping"]


def family_difference_text(family_id: str) -> str:
    return FAMILY_DEFAULTS.get(family_id, FAMILY_DEFAULTS["util_full_stack"])["difference"]


def family_next_action_text(family_id: str, status: str) -> str:
    base = FAMILY_DEFAULTS.get(family_id, FAMILY_DEFAULTS["util_full_stack"])["next_action"]
    if status == STATUS_EQUIVALENT:
        return "保持现有代码路径、测试和 gate 证据不回退。"
    if status == STATUS_AGGREGATED:
        return base
    if status == STATUS_APPROX:
        return base
    if status == STATUS_MISSING:
        return "需要补显式代码落点、config/route 入口与至少一个 smoke 或 matrix case。"
    return "保持为非主线对齐对象；无需纳入当前 Surface Charging 审计主线。"


def selectors_to_strings(test_evidence: list[dict]) -> list[str]:
    flattened: list[str] = []
    for item in test_evidence:
        path = norm_path_string(item["path"])
        selectors = item.get("selectors", [])
        if selectors:
            flattened.append(f"{path}::{', '.join(selectors)}")
        else:
            flattened.append(path)
    return flattened


def main() -> None:
    family_matrix = load_json(FAMILY_MATRIX_PATH)
    sealoff = load_json(SEALOFF_PATH)

    target_families = {item["id"]: item for item in family_matrix["target_families"]}
    family_to_code = {item["family_id"]: item for item in family_matrix["family_to_code"]}
    family_to_test = {item["family_id"]: item for item in family_matrix["family_to_test"]}
    family_to_config = {item["family_id"]: item for item in family_matrix["family_to_config_route"]}

    entries = []
    for java_path in sorted(SPIS_ROOT.rglob("*.java")):
        relative = java_path.relative_to(SPIS_ROOT).as_posix()
        family_id, family_name = resolve_family(relative)
        family_target = target_families.get(family_id, {"domain": relative.split("/")[0], "spis_family": family_name})
        status = infer_status(relative, family_id)

        code_info = family_to_code.get(family_id, {"code_paths": [], "documentation_paths": []})
        test_info = family_to_test.get(family_id, {"test_evidence": [], "gate_artifacts": []})
        config_info = family_to_config.get(
            family_id,
            {
                "metadata_keys": [],
                "metadata_source_files": [],
                "config_keys": [],
                "config_source_files": [],
                "preset_names": [],
                "reference_matrix_cases": [],
            },
        )

        notes = []
        if status in {STATUS_AGGREGATED, STATUS_APPROX} and family_target.get("alignment_status") == "completed":
            notes.append("family 已封板；当前结论属于 class 级聚合承载差异，不直接推翻 sealoff。")
        if family_id in {"util_func", "util_instrument", "top_default"}:
            notes.append("该 family 在当前项目中以公共基础设施或运行时对象层为主，类级镜像天然偏弱。")

        entry = {
            "spis_class_path": relative,
            "spis_domain": family_target["domain"],
            "spis_family": family_name,
            "family_id": family_id,
            "current_project_mapping": family_mapping_text(family_id, relative),
            "mapping_status": status,
            "evidence": {
                "code_paths": code_info.get("code_paths", []),
                "config_route_paths": sorted(
                    {norm_path_string(item) for item in config_info.get("metadata_source_files", [])}
                    | {norm_path_string(item) for item in config_info.get("config_source_files", [])}
                    | {rel(REFERENCE_MATRIX_PATH)}
                ),
                "test_paths": selectors_to_strings(test_info.get("test_evidence", [])),
                "gate_artifacts": sorted(
                    {norm_path_string(item) for item in test_info.get("gate_artifacts", [])}
                    | {
                        rel(GATE_PATH),
                        rel(SEALOFF_PATH),
                    }
                ),
                "documentation_paths": sorted(
                    {norm_path_string(item) for item in code_info.get("documentation_paths", [])}
                    | {
                        rel(PLAN_DOC_PATH),
                        rel(CONTRACT_PATH),
                    }
                ),
            },
            "key_difference": family_difference_text(family_id),
            "next_action": family_next_action_text(family_id, status),
            "notes": notes,
        }
        entries.append(entry)

    domain_status_counts: dict[str, Counter] = defaultdict(Counter)
    family_status_counts: dict[str, Counter] = defaultdict(Counter)
    family_counts: dict[str, int] = defaultdict(int)
    for entry in entries:
        domain_status_counts[entry["spis_domain"]][entry["mapping_status"]] += 1
        family_status_counts[entry["family_id"]][entry["mapping_status"]] += 1
        family_counts[entry["family_id"]] += 1

    backlog = []
    for family_id, counts in sorted(family_status_counts.items()):
        approx_or_missing = counts[STATUS_APPROX] + counts[STATUS_MISSING]
        if approx_or_missing == 0:
            continue
        backlog.append(
            {
                "family_id": family_id,
                "spis_family": target_families.get(family_id, {}).get("spis_family", family_id),
                "approximate_or_missing_class_count": approx_or_missing,
                "recommendation": FAMILY_ATTENTION_NOTES.get(
                    family_id,
                    "优先补类级职责索引、显式对象层映射和更直接的单测/route 证据。",
                ),
            }
        )

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "scope": "Surface Charging 主线边界内六大域逐 class 审计",
        "source_root": str(SPIS_ROOT.relative_to(REPO_ROOT)),
        "tracked_family_count": sealoff["tracked_family_count"],
        "sealoff_alignment_status_counts": sealoff["alignment_status_counts"],
        "final_alignment_verdict": sealoff["final_alignment_verdict"],
        "total_spis_class_count": len(entries),
        "domain_status_counts": {domain: dict(counter) for domain, counter in sorted(domain_status_counts.items())},
        "family_status_counts": {family: dict(counter) for family, counter in sorted(family_status_counts.items())},
        "backlog": backlog,
    }

    CLASS_AUDIT_JSON.write_text(
        json.dumps(
            {
                "metadata": {
                    "generated_at_utc": summary["generated_at_utc"],
                    "scope": summary["scope"],
                    "source_root": summary["source_root"],
                    "baseline_artifacts": {
                        "plan_doc": rel(PLAN_DOC_PATH),
                        "family_matrix": rel(FAMILY_MATRIX_PATH),
                        "contract": rel(CONTRACT_PATH),
                        "coverage_gate": rel(GATE_PATH),
                        "sealoff": rel(SEALOFF_PATH),
                        "reference_matrix": rel(REFERENCE_MATRIX_PATH),
                    },
                    "sealoff_baseline": {
                        "tracked_family_count": sealoff["tracked_family_count"],
                        "alignment_status_counts": sealoff["alignment_status_counts"],
                        "final_alignment_verdict": sealoff["final_alignment_verdict"],
                    },
                },
                "summary": summary,
                "entries": entries,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    with CLASS_AUDIT_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "spis_class_path",
                "spis_domain",
                "spis_family",
                "family_id",
                "current_project_mapping",
                "mapping_status",
                "evidence_code_paths",
                "evidence_config_route_paths",
                "evidence_test_paths",
                "evidence_gate_artifacts",
                "evidence_documentation_paths",
                "key_difference",
                "next_action",
                "notes",
            ]
        )
        for entry in entries:
            writer.writerow(
                [
                    entry["spis_class_path"],
                    entry["spis_domain"],
                    entry["spis_family"],
                    entry["family_id"],
                    entry["current_project_mapping"],
                    entry["mapping_status"],
                    " | ".join(entry["evidence"]["code_paths"]),
                    " | ".join(entry["evidence"]["config_route_paths"]),
                    " | ".join(entry["evidence"]["test_paths"]),
                    " | ".join(entry["evidence"]["gate_artifacts"]),
                    " | ".join(entry["evidence"]["documentation_paths"]),
                    entry["key_difference"],
                    entry["next_action"],
                    " | ".join(entry["notes"]),
                ]
            )

    domain_table_lines = []
    ordered_statuses = [STATUS_EQUIVALENT, STATUS_AGGREGATED, STATUS_APPROX, STATUS_MISSING, STATUS_OUT]
    for domain in ["Surf", "Vol", "Top", "Circ", "Solver", "Util", "osgi"]:
        if domain not in domain_status_counts:
            continue
        counter = domain_status_counts[domain]
        total = sum(counter.values())
        domain_table_lines.append(
            "| {domain} | {total} | {equiv} | {agg} | {approx} | {missing} | {out} |".format(
                domain=domain,
                total=total,
                equiv=counter.get(STATUS_EQUIVALENT, 0),
                agg=counter.get(STATUS_AGGREGATED, 0),
                approx=counter.get(STATUS_APPROX, 0),
                missing=counter.get(STATUS_MISSING, 0),
                out=counter.get(STATUS_OUT, 0),
            )
        )

    backlog_lines = []
    for item in backlog[:12]:
        backlog_lines.append(
            f"- `{item['family_id']}`: `{item['approximate_or_missing_class_count']}` 个 class 仍是 `近似实现/尚无实现`。{item['recommendation']}"
        )
    if not backlog_lines:
        backlog_lines.append("- 当前没有 `近似实现` 或 `尚无实现` 项。")

    md = f"""# Surface/SPIS 逐 Class 审计摘要

## 1. 基线

- 生成时间：`{summary["generated_at_utc"]}`
- 审计范围：`{summary["scope"]}`
- SPIS class 总数：`{summary["total_spis_class_count"]}`
- sealoff 基线：`tracked_family_count = {sealoff["tracked_family_count"]}`，`final_alignment_verdict = {sealoff["final_alignment_verdict"]}`
- 参考主文档：`{rel(PLAN_DOC_PATH)}`

## 2. 产物

- 机器可读主表：`{rel(CLASS_AUDIT_JSON)}`
- 可筛选 CSV：`{rel(CLASS_AUDIT_CSV)}`
- family matrix：`{rel(FAMILY_MATRIX_PATH)}`
- coverage gate：`{rel(GATE_PATH)}`
- sealoff：`{rel(SEALOFF_PATH)}`

## 3. 域级统计

| 域 | class 数 | 已等价实现 | 聚合等价实现 | 近似实现 | 尚无实现 | 不纳入范围 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{chr(10).join(domain_table_lines)}

## 4. 解读规则

- `已等价实现`：当前项目里有较明确的单类或强命名落点，且证据链可回溯。
- `聚合等价实现`：语义已闭环，但由多个 C++ 模块共同承载，不是 Java 单类镜像。
- `近似实现`：family 已有 route/metadata/gate 证据，但 class 级对象层或独立语义承载偏弱。
- `尚无实现`：当前边界内未找到足够证据。
- `不纳入当前目标范围`：当前审计不纳入 Surface Charging 主线。

> 说明：如果某个 class 被标为 `聚合等价实现` 或 `近似实现`，并不自动推翻当前 family sealoff；这表示“family 已封板，但 class 级对象层仍有聚合承载差异”。

## 5. 后续任务摘要

{chr(10).join(backlog_lines)}

## 6. 高风险热点

- `Toolkit/Surface Charging/src/SurfaceChargingKernelFramework.cpp`
- `Toolkit/Surface Charging/src/DensePlasmaSurfaceCharging.cpp`
- `Tools/FieldSolver/src/SurfaceFieldVolumeBridge.cpp`

这些文件仍是 class-level 语义被集中承载的主要热点，后续如果继续追求更强的逐 class 可审计性，应优先从这里拆对象层索引。
"""
    CLASS_AUDIT_SUMMARY.write_text(md, encoding="utf-8")


if __name__ == "__main__":
    main()
