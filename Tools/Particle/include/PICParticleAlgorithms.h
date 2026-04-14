#ifndef SCDAT_PIC_PARTICLE_ALGORITHMS_H
#define SCDAT_PIC_PARTICLE_ALGORITHMS_H

#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"
#include "../../Geometry/include/Matrix3x3.h"
#include "../../Mesh/include/MeshAlgorithms.h"
#include "../../Mesh/include/MeshPartitioning.h"
#include "../../Mesh/include/MeshParsing.h"
#include "ParticleTypes.h"
#include "ParticleManager.h"
#include <atomic>
#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

namespace SCDAT {
namespace Particle {

// 鍓嶅悜澹版槑
class CrossTetrahedron;
class FieldInterpolator;
class ChargeDepositor;
using CrossTetrahedronPtr = std::shared_ptr<CrossTetrahedron>;
using FieldInterpolatorPtr = std::shared_ptr<FieldInterpolator>;
using ChargeDepositorPtr = std::shared_ptr<ChargeDepositor>;

using NodeId = Mesh::NodeId;
using ElementId = Mesh::ElementId;

enum class TrajectoryChargeDepositionKernel
{
    NearestPoint,
    LinearSegmentCloud,
    QuadratureGauss3
};

/**
 * @brief 绮掑瓙浣嶇疆璺熻釜缁撴瀯 / Particle Location Tracking Structure
 * @details 鐢ㄤ簬璺熻釜绮掑瓙鍦ㄩ潪缁撴瀯鍖栫綉鏍间腑鐨勭簿纭綅缃俊鎭紝
 *          鍖呮嫭鎵€鍦ㄥ崟鍏冨拰閲嶅績鍧愭爣
 *
 * 鏁版嵁缁撴瀯 / Data Structure:
 * - 绮掑瓙鎸囬拡锛氭寚鍚戣璺熻釜鐨勭矑瀛愬璞?
 * - 褰撳墠鍗曞厓锛氱矑瀛愭墍鍦ㄧ殑缃戞牸鍗曞厓ID
 * - 閲嶅績鍧愭爣锛氱矑瀛愬湪鍗曞厓鍐呯殑閲嶅績鍧愭爣
 * - 鏈夋晥鎬ф爣蹇楋細浣嶇疆淇℃伅鏄惁鏈夋晥
 *
 * 搴旂敤鍦烘櫙 / Use Cases:
 * - 绮掑瓙杞ㄨ抗绉垎涓殑浣嶇疆璺熻釜
 * - 鍦烘彃鍊兼椂鐨勫揩閫熷畾浣?
 * - 杈圭晫鏉′欢澶勭悊
 * - 绮掑瓙涓㈠け妫€娴?
 *
 * @note 閲嶅績鍧愭爣婊¤冻锛氣垜位岬?= 1 涓?位岬?鈮?0
 * @warning 褰撶矑瀛愮Щ鍑鸿绠楀煙鏃讹紝is_valid灏嗚璁剧疆涓篺alse
 */
struct ParticleLocation {
    ParticlePtr particle;                    ///< 绮掑瓙鎸囬拡
    ElementId current_element;               ///< 褰撳墠鎵€鍦ㄥ崟鍏僆D
    std::vector<double> barycentric_coords;  ///< 閲嶅績鍧愭爣锛堝洓闈綋锛?涓垎閲忥級
    bool is_valid;                          ///< 浣嶇疆淇℃伅鏈夋晥鎬ф爣蹇?

    /**
     * @brief 榛樿鏋勯€犲嚱鏁?
     * @details 鍒涘缓鏃犳晥鐨勭矑瀛愪綅缃璞?
     */
    ParticleLocation() : current_element(0), is_valid(false) {}

    /**
     * @brief 鏋勯€犲嚱鏁?
     * @param p 绮掑瓙鎸囬拡
     * @param elem 鍗曞厓ID
     * @details 鍒涘缓鏈夋晥鐨勭矑瀛愪綅缃璞★紝閲嶅績鍧愭爣闇€瑕佸悗缁绠?
     */
    ParticleLocation(ParticlePtr p, ElementId elem)
        : particle(p), current_element(elem), is_valid(true) {}
};

/**
 * @brief 鍦烘彃鍊煎櫒绫?
 *
 * 璐熻矗灏嗙綉鏍间笂鐨勫満鍊兼彃鍊煎埌绮掑瓙浣嶇疆
 * 鍩轰簬SPIS鍘熺増鐨勫満鎻掑€肩畻娉?
 */
class FieldInterpolator {
public:
    // 鏋勯€犲嚱鏁?
    FieldInterpolator();
    ~FieldInterpolator() = default;

    // 璁剧疆缃戞牸
    void setMesh(const std::vector<Mesh::NodePtr>& nodes, const std::vector<Mesh::ElementPtr>& elements);
    void setNodes(const std::vector<Mesh::NodePtr>& nodes) { nodes_ = nodes; }
    void setElements(const std::vector<Mesh::ElementPtr>& elements) { elements_ = elements; }

    // 璁剧疆鍦烘暟鎹?
    void setElectricField(const std::vector<Utils::Vector3D>& electric_field);
    void setMagneticField(const std::vector<Utils::Vector3D>& magnetic_field);
    void setPotential(const std::vector<double>& potential);

    // 鍦烘彃鍊兼柟娉?
    Utils::Vector3D interpolateElectricField(const Utils::Point3D& position, ElementId element_id) const;
    Utils::Vector3D interpolateMagneticField(const Utils::Point3D& position, ElementId element_id) const;
    double interpolatePotential(const Utils::Point3D& position, ElementId element_id) const;

    // 姊害璁＄畻
    Utils::Vector3D computeElectricFieldGradient(const Utils::Point3D& position, ElementId element_id) const;

    // 閲嶅績鍧愭爣璁＄畻
    bool computeBarycentricCoordinates(const Utils::Point3D& position, ElementId element_id,
                                      std::vector<double>& barycentric) const;

private:
    std::vector<Mesh::NodePtr> nodes_;
    std::vector<Mesh::ElementPtr> elements_;
    std::vector<Utils::Vector3D> electric_field_;
    std::vector<Utils::Vector3D> magnetic_field_;
    std::vector<double> potential_;

    // 杈呭姪鏂规硶
    bool isPointInTetrahedron(const Utils::Point3D& point, ElementId element_id,
                             std::vector<double>& barycentric) const;
    Utils::Vector3D interpolateVectorField(const Utils::Point3D& position, ElementId element_id,
                                          const std::vector<Utils::Vector3D>& field) const;
    double interpolateScalarField(const Utils::Point3D& position, ElementId element_id,
                                 const std::vector<double>& field) const;

    // 鍑犱綍棰勮绠楃紦瀛?
    struct ElementGeometry {
        Utils::Point3D origin;        // P0
        Utils::Matrix3x3 inv_jacobian; // J^-1
        bool is_valid = false;
    };
    std::vector<ElementGeometry> geometry_cache_;
    void precomputeGeometry();
};

/**
 * @brief 鐢佃嵎娌夌Н鍣ㄧ被
 *
 * 璐熻矗灏嗙矑瀛愮數鑽锋矇绉埌缃戞牸鑺傜偣
 * 瀹炵幇鐢佃嵎瀹堟亽鐨勬矇绉畻娉?
 */
class ChargeDepositor {
public:
    // 鏋勯€犲嚱鏁?
    ChargeDepositor();
    ~ChargeDepositor() = default;

    // 璁剧疆缃戞牸
    void setMesh(const std::vector<Mesh::NodePtr>& nodes, const std::vector<Mesh::ElementPtr>& elements);

    // 鐢佃嵎娌夌Н鏂规硶
    void depositCharge(const Utils::Point3D& position, ElementId element_id,
                      double charge, double weight, std::vector<double>& charge_density);

    // 杞ㄨ抗绉垎鐢佃嵎娌夌Н
    void depositChargeAlongTrajectory(const Utils::Point3D& start_pos, const Utils::Point3D& end_pos,
                                     ElementId element_id, double charge, double weight, double dt,
                                     std::vector<double>& charge_density);

    void setTrajectoryDepositionScheme(TrajectoryChargeDepositionKernel kernel,
                                       std::size_t segment_count);

    // 娓呴浂鐢佃嵎瀵嗗害
    void clearChargeDensity(std::vector<double>& charge_density);

    // 鐢佃嵎瀹堟亽妫€鏌?
    double computeTotalCharge(const std::vector<double>& charge_density) const;

private:
    std::vector<Mesh::NodePtr> nodes_;
    std::vector<Mesh::ElementPtr> elements_;
    std::shared_ptr<FieldInterpolator> interpolator_; // Cached interpolator

    // 杈呭姪鏂规硶
    void performDeposit(const std::vector<double>& barycentric, const std::vector<NodeId>& nodes,
                       double charge_contribution, std::vector<double>& charge_density);

    TrajectoryChargeDepositionKernel trajectory_kernel_ =
        TrajectoryChargeDepositionKernel::LinearSegmentCloud;
    std::size_t trajectory_segment_count_ = 10;
};

/**
 * @brief 鍥涢潰浣撶┛瓒婄畻娉曠被
 *
 * 鍩轰簬SPIS鍘熺増CrossTetrahedron绠楁硶鐨凜++瀹炵幇
 * 澶勭悊绮掑瓙鍦ㄥ洓闈綋缃戞牸涓殑杞ㄨ抗绉垎
 */
class CrossTetrahedron {
public:
    // 杞ㄨ抗绉垎缁撴灉
    struct TrajectoryResult {
        Utils::Point3D final_position;     ///< 鏈€缁堜綅缃?
        Utils::Vector3D final_velocity;    ///< 鏈€缁堥€熷害
        ElementId final_element;           ///< 鏈€缁堟墍鍦ㄥ崟鍏?
        double integration_time;           ///< 瀹為檯绉垎鏃堕棿
        bool crossed_boundary;             ///< 鏄惁绌胯秺杈圭晫
        bool lost_particle;                ///< 鏄惁涓㈠け绮掑瓙
        std::string loss_reason;           ///< 涓㈠け鍘熷洜
    };

    // 鏋勯€犲嚱鏁?
    CrossTetrahedron();
    ~CrossTetrahedron() = default;

    // 璁剧疆缃戞牸鍜屽満
    void setMesh(const std::vector<Mesh::NodePtr>& nodes, const std::vector<Mesh::ElementPtr>& elements);
    void setFieldInterpolator(FieldInterpolatorPtr interpolator);
    void setChargeDepositor(ChargeDepositorPtr depositor);

    // 涓昏杞ㄨ抗绉垎鏂规硶
    TrajectoryResult crossTetrahedron(ParticlePtr particle, double dt,
                                     std::vector<double>& charge_density,
                                     bool integrate_deposit = true,
                                     ElementId element_hint = static_cast<ElementId>(-1));

    // 涓嶅悓鍦洪厤缃殑杞ㄨ抗绉垎
    TrajectoryResult crossTetraUniformE(ParticlePtr particle, double dt,
                                       const Utils::Vector3D& uniform_E,
                                       std::vector<double>& charge_density);

    TrajectoryResult crossTetraVaryingE(ParticlePtr particle, double dt,
                                       std::vector<double>& charge_density);

    TrajectoryResult crossTetraUniformEwithB(ParticlePtr particle, double dt,
                                            const Utils::Vector3D& uniform_E,
                                            const Utils::Vector3D& uniform_B,
                                            std::vector<double>& charge_density);

    // 杈圭晫澶勭悊
    void setBoundaryReflection(bool enable, double reflection_potential = 0.0);
    bool handleBoundaryCondition(ParticlePtr particle, const Utils::Point3D& intersection_point);

    // 缁熻淇℃伅
    struct Statistics {
        int total_crossings = 0;
        int boundary_crossings = 0;
        int lost_particles = 0;
        int reflection_events = 0;
        double total_integration_time = 0.0;
    };

    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics() { stats_ = Statistics(); }

private:
    std::vector<Mesh::NodePtr> nodes_;
    std::vector<Mesh::ElementPtr> elements_;
    FieldInterpolatorPtr field_interpolator_;
    ChargeDepositorPtr charge_depositor_;

    // 杈圭晫鏉′欢鍙傛暟
    bool boundary_reflection_;
    double reflection_potential_;

    // 缁熻淇℃伅
    Statistics stats_;

    // 杞ㄨ抗绉垎鏍稿績绠楁硶
    TrajectoryResult integrateTrajectory(ParticlePtr particle, double dt,
                                        std::vector<double>& charge_density,
                                        bool uniform_field = false,
                                        const Utils::Vector3D& E_uniform = Utils::Vector3D(),
                                        const Utils::Vector3D& B_uniform = Utils::Vector3D(),
                                        ElementId element_hint = static_cast<ElementId>(-1));

    // 杈圭晫妫€娴?
    bool detectBoundaryCrossing(const Utils::Point3D& start, const Utils::Point3D& end,
                               ElementId element_id, Utils::Point3D& intersection);

    // 鏃堕棿姝ラ暱鎺у埗
    double adaptiveTimeStep(ParticlePtr particle, double dt_max, const Utils::Vector3D& E, const Utils::Vector3D& B);

    // Boris绠楁硶绮掑瓙鎺ㄨ繘
    Utils::Vector3D borisIntegration(const Utils::Vector3D& velocity,
                                    const Utils::Vector3D& E, const Utils::Vector3D& B,
                                    double charge_over_mass, double dt);

    // 绮掑瓙瀹氫綅
    ElementId locateParticle(const Utils::Point3D& position, ElementId hint = static_cast<ElementId>(-1));

    // 浜屽垎娉曟眰瑙ｈ竟鐣屼氦鐐?
    bool dichotomyBoundaryIntersection(const Utils::Point3D& start, const Utils::Point3D& end,
                                      ElementId element_id, Utils::Point3D& intersection, double tolerance = 1e-12);
};

} // namespace Particle
} // namespace SCDAT

#endif // SCDAT_PIC_PARTICLE_ALGORITHMS_H
