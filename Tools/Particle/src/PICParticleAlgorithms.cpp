#include "PICParticleAlgorithms.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace SCDAT
{
namespace Particle
{

//==============================================================================
// FieldInterpolator 瀹炵幇
//==============================================================================

FieldInterpolator::FieldInterpolator() {}

void FieldInterpolator::setMesh(const std::vector<Mesh::NodePtr>& nodes,
                                const std::vector<Mesh::ElementPtr>& elements)
{
    nodes_ = nodes;
    elements_ = elements;
    precomputeGeometry();
}

void FieldInterpolator::precomputeGeometry()
{
    if (elements_.empty())
        return;

    geometry_cache_.resize(elements_.size());

// IMPORTANT: 绂佺敤 Jacobian 閫嗘柟娉曪紝寮哄埗浣跨敤 signed volume 鏂规硶
// 鍘熷洜锛欽acobian 閫嗘柟娉曞湪鏌愪簺缃戞牸閰嶇疆涓嬩笉鍑嗙‘锛屼笌 TetrahedronElement::contains() 淇淇濇寔涓€锟?
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < static_cast<int>(elements_.size()); ++i)
    {
        geometry_cache_[i].is_valid = false; // 绂佺敤缂撳瓨锛屽己鍒朵娇锟?fallback volume method
    }

    // 鍘熷 Jacobian 閫嗘柟娉曚唬鐮佸凡琚敞閲婏紙瀛樺湪鏁板€间笉绋冲畾闂锟?
    /*
    for (int i = 0; i < static_cast<int>(elements_.size()); ++i) {
        const auto& element = elements_[i];
        const auto& nodes = element->getNodes();

        if (nodes.size() != 4) {
            geometry_cache_[i].is_valid = false;
            continue;
        }

        const auto& p0 = nodes[0]->getPosition();
        const auto& p1 = nodes[1]->getPosition();
        const auto& p2 = nodes[2]->getPosition();
        const auto& p3 = nodes[3]->getPosition();

        Utils::Vector3D v1 = p1 - p0;
        Utils::Vector3D v2 = p2 - p0;
        Utils::Vector3D v3 = p3 - p0;

        // J = [v1, v2, v3]
        Utils::Matrix3x3 J = Utils::Matrix3x3::fromColumns(v1, v2, v3);

        if (std::abs(J.determinant()) > 1e-15) {
            geometry_cache_[i].inv_jacobian = J.inverse();
            geometry_cache_[i].origin = p0;
            geometry_cache_[i].is_valid = true;
        } else {
            geometry_cache_[i].is_valid = false;
        }
    }
    */
}

void FieldInterpolator::setElectricField(const std::vector<Utils::Vector3D>& electric_field)
{
    electric_field_ = electric_field;
}

void FieldInterpolator::setMagneticField(const std::vector<Utils::Vector3D>& magnetic_field)
{
    magnetic_field_ = magnetic_field;
}

void FieldInterpolator::setPotential(const std::vector<double>& potential)
{
    potential_ = potential;
}

Utils::Vector3D FieldInterpolator::interpolateElectricField(const Utils::Point3D& position,
                                                            ElementId element_id) const
{
    if (elements_.empty() || electric_field_.empty())
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    return interpolateVectorField(position, element_id, electric_field_);
}

Utils::Vector3D FieldInterpolator::interpolateMagneticField(const Utils::Point3D& position,
                                                            ElementId element_id) const
{
    if (elements_.empty() || magnetic_field_.empty())
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    return interpolateVectorField(position, element_id, magnetic_field_);
}

double FieldInterpolator::interpolatePotential(const Utils::Point3D& position,
                                               ElementId element_id) const
{
    if (elements_.empty() || potential_.empty())
    {
        return 0.0;
    }

    return interpolateScalarField(position, element_id, potential_);
}

bool FieldInterpolator::computeBarycentricCoordinates(const Utils::Point3D& position,
                                                      ElementId element_id,
                                                      std::vector<double>& barycentric) const
{
    if (elements_.empty())
        return false;

    return isPointInTetrahedron(position, element_id, barycentric);
}

Utils::Vector3D
FieldInterpolator::interpolateVectorField(const Utils::Point3D& position, ElementId element_id,
                                          const std::vector<Utils::Vector3D>& field) const
{
    std::vector<double> barycentric(4);
    if (!computeBarycentricCoordinates(position, element_id, barycentric))
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    // 鑾峰彇鍥涢潰浣撶殑鑺傜偣
    if (element_id >= elements_.size())
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    const auto& element = elements_[element_id];
    const auto& nodes = element->getNodes();

    // 绾挎€ф彃锟?
    Utils::Vector3D result(0.0, 0.0, 0.0);
    for (size_t i = 0; i < 4 && i < nodes.size(); ++i)
    {
        NodeId node_id = nodes[i]->getId();
        if (node_id < field.size())
        {
            result = result + field[node_id] * barycentric[i];
        }
    }

    return result;
}

double FieldInterpolator::interpolateScalarField(const Utils::Point3D& position,
                                                 ElementId element_id,
                                                 const std::vector<double>& field) const
{
    std::vector<double> barycentric(4);
    if (!computeBarycentricCoordinates(position, element_id, barycentric))
    {
        return 0.0;
    }

    // 鑾峰彇鍥涢潰浣撶殑鑺傜偣
    if (element_id >= elements_.size())
    {
        return 0.0;
    }

    const auto& element = elements_[element_id];
    const auto& nodes = element->getNodes();

    // 绾挎€ф彃锟?
    double result = 0.0;
    for (size_t i = 0; i < 4 && i < nodes.size(); ++i)
    {
        NodeId node_id = nodes[i]->getId();
        if (node_id < field.size())
        {
            result += field[node_id] * barycentric[i];
        }
    }

    return result;
}

bool FieldInterpolator::isPointInTetrahedron(const Utils::Point3D& point, ElementId element_id,
                                             std::vector<double>& barycentric) const
{
    if (elements_.empty() || element_id >= elements_.size())
        return false;

    // Optimization: Use precomputed inverse Jacobian
    if (element_id < geometry_cache_.size() && geometry_cache_[element_id].is_valid)
    {
        const auto& geom = geometry_cache_[element_id];

        // P - P0
        Utils::Vector3D dP = point - geom.origin;

        // lambda = J^-1 * (P - P0)
        Utils::Vector3D lambda = geom.inv_jacobian * dP;

        barycentric[1] = lambda.x();
        barycentric[2] = lambda.y();
        barycentric[3] = lambda.z();
        barycentric[0] = 1.0 - lambda.x() - lambda.y() - lambda.z();

        // FIX: Correct epsilon logic with proper bounds checking
        // For 1D degenerate tetrahedra, numerical errors require larger tolerance
        const double epsilon = 1e-8; // Increased from 1e-10 for robustness
        return (barycentric[0] >= -epsilon && barycentric[0] <= 1.0 + epsilon &&
                barycentric[1] >= -epsilon && barycentric[1] <= 1.0 + epsilon &&
                barycentric[2] >= -epsilon && barycentric[2] <= 1.0 + epsilon &&
                barycentric[3] >= -epsilon && barycentric[3] <= 1.0 + epsilon);
    }

    // Fallback to volume method if cache is invalid (e.g. degenerate element)
    const auto& element = elements_[element_id];
    const auto& nodes = element->getNodes();

    // 鑾峰彇鍥涢潰浣撶殑鍥涗釜椤剁偣
    // Optimization: Use std::array instead of std::vector
    std::array<Utils::Point3D, 4> vertices;
    for (size_t i = 0; i < 4 && i < nodes.size(); ++i)
    {
        vertices[i] = nodes[i]->getPosition();
    }

    // 璁＄畻閲嶅績鍧愭爣 - 浣跨敤 SIGNED volume锛堜笉浣跨敤 abs()锟?
    // 杩欎笌 TetrahedronElement::contains() 涓慨澶嶅悗鐨勭畻娉曚竴锟?
    auto signedVolume6 = [](const Utils::Point3D& a, const Utils::Point3D& b,
                            const Utils::Point3D& c, const Utils::Point3D& d) -> double
    {
        Utils::Vector3D v1 = b - a;
        Utils::Vector3D v2 = c - a;
        Utils::Vector3D v3 = d - a;
        return v1.dot(v2.cross(v3)); // 鏍囬噺涓夐噸绉紝淇濈暀绗﹀彿
    };

    double V6_total = signedVolume6(vertices[0], vertices[1], vertices[2], vertices[3]);

    if (std::abs(V6_total) < 1e-15)
    {
        // 閫€鍖栧洓闈綋
        return false;
    }

    // 璁＄畻鍚勪釜瀛愬洓闈綋锟?signed volume锛岄櫎浠ワ拷?volume 寰楀埌閲嶅績鍧愭爣
    barycentric[0] = signedVolume6(point, vertices[1], vertices[2], vertices[3]) / V6_total;
    barycentric[1] = signedVolume6(vertices[0], point, vertices[2], vertices[3]) / V6_total;
    barycentric[2] = signedVolume6(vertices[0], vertices[1], point, vertices[3]) / V6_total;
    barycentric[3] = signedVolume6(vertices[0], vertices[1], vertices[2], point) / V6_total;

    // 妫€鏌ユ槸鍚﹀湪鍥涢潰浣撳唴閮紙鎵€鏈夐噸蹇冨潗鏍囬兘锟?[-epsilon, 1+epsilon] 鑼冨洿鍐咃級
    // FIX: Increased epsilon for 1D degenerate tetrahedra numerical stability
    const double epsilon = 1e-8; // Increased from 1e-10 to handle numerical errors
    for (double coord : barycentric)
    {
        if (coord < -epsilon || coord > 1.0 + epsilon)
        {
            return false;
        }
    }

    return true;
}

//==============================================================================
// ChargeDepositor 瀹炵幇
//==============================================================================

ChargeDepositor::ChargeDepositor()
{
    interpolator_ = std::make_shared<FieldInterpolator>();
}

void ChargeDepositor::setMesh(const std::vector<Mesh::NodePtr>& nodes,
                              const std::vector<Mesh::ElementPtr>& elements)
{
    nodes_ = nodes;
    elements_ = elements;
    if (interpolator_)
    {
        interpolator_->setMesh(nodes, elements);
    }
}

void ChargeDepositor::setTrajectoryDepositionScheme(
    TrajectoryChargeDepositionKernel kernel, std::size_t segment_count)
{
    trajectory_kernel_ = kernel;
    trajectory_segment_count_ = std::clamp<std::size_t>(segment_count, 1, 128);
}

void ChargeDepositor::depositCharge(const Utils::Point3D& position, ElementId element_id,
                                    double charge, double weight,
                                    std::vector<double>& charge_density)
{
    if (elements_.empty() || !interpolator_)
        return;

    // 璁＄畻閲嶅績鍧愭爣
    // Optimization: Use std::array or hoisted vector to avoid allocation
    // Since computeBarycentricCoordinates takes vector&, we use a thread_local or static if
    // possible, but for now let's just hoist it if we can, or accept one small allocation vs the
    // huge copy of mesh before. Actually, let's use a local vector but at least we saved the mesh
    // copy.
    std::vector<double> barycentric(4);

    if (!interpolator_->computeBarycentricCoordinates(position, element_id, barycentric))
    {
        return; // 鐐逛笉鍦ㄥ洓闈綋锟?
    }

    // 鑾峰彇鍥涢潰浣撶殑鑺傜偣
    if (element_id >= elements_.size())
        return;
    const auto& element = elements_[element_id];
    const auto& nodes = element->getNodes();

    // 鎵ц鐢佃嵎娌夌Н
    // Optimization: Avoid vector allocation for node_ids
    std::array<NodeId, 4> node_ids;
    size_t count = 0;
    for (const auto& node : nodes)
    {
        if (count < 4)
            node_ids[count++] = node->getId();
    }

    // Inline performDeposit logic to avoid vector creation
    for (size_t i = 0; i < count; ++i)
    {
        if (node_ids[i] < charge_density.size())
        {
            charge_density[node_ids[i]] += (charge * weight) * barycentric[i];
        }
    }
}

void ChargeDepositor::depositChargeAlongTrajectory(const Utils::Point3D& start_pos,
                                                   const Utils::Point3D& end_pos,
                                                   ElementId element_id, double charge,
                                                   double weight, double /* dt */,
                                                   std::vector<double>& charge_density)
{
    if (elements_.empty())
        return;

    Utils::Vector3D trajectory = end_pos - start_pos;
    if (trajectory.magnitude() < 1e-15)
    {
        depositCharge(start_pos, element_id, charge, weight, charge_density);
        return;
    }

    if (trajectory_kernel_ == TrajectoryChargeDepositionKernel::NearestPoint)
    {
        depositCharge(end_pos, element_id, charge, weight, charge_density);
        return;
    }

    const std::size_t num_segments = std::max<std::size_t>(1, trajectory_segment_count_);

    if (trajectory_kernel_ == TrajectoryChargeDepositionKernel::QuadratureGauss3)
    {
        // Gauss-Legendre 3-point quadrature on each segment for smoother deposition.
        constexpr std::array<double, 3> kNodes = {
            -0.7745966692414834, 0.0, 0.7745966692414834};
        constexpr std::array<double, 3> kWeights = {
            5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};

        for (std::size_t segment = 0; segment < num_segments; ++segment)
        {
            const double segment_start_t = static_cast<double>(segment) /
                                           static_cast<double>(num_segments);
            const double segment_end_t = static_cast<double>(segment + 1) /
                                         static_cast<double>(num_segments);
            const double segment_center_t = 0.5 * (segment_start_t + segment_end_t);
            const double segment_half_span_t = 0.5 * (segment_end_t - segment_start_t);

            for (std::size_t q = 0; q < kNodes.size(); ++q)
            {
                const double t = segment_center_t + segment_half_span_t * kNodes[q];
                const Utils::Point3D sample_point = start_pos + trajectory * t;
                const double weighted_charge = charge * kWeights[q] /
                                               static_cast<double>(num_segments);
                depositCharge(sample_point, element_id, weighted_charge, weight, charge_density);
            }
        }
        return;
    }

    // Linear cloud-in-cell style segment-center deposition.
    const double charge_per_segment = charge / static_cast<double>(num_segments);
    for (std::size_t segment = 0; segment < num_segments; ++segment)
    {
        const double t = (static_cast<double>(segment) + 0.5) /
                         static_cast<double>(num_segments);
        const Utils::Point3D segment_center = start_pos + trajectory * t;
        depositCharge(segment_center, element_id, charge_per_segment, weight, charge_density);
    }
}

void ChargeDepositor::clearChargeDensity(std::vector<double>& charge_density)
{
    std::fill(charge_density.begin(), charge_density.end(), 0.0);
}

double ChargeDepositor::computeTotalCharge(const std::vector<double>& charge_density) const
{
    double total = 0.0;
    for (double charge : charge_density)
    {
        total += charge;
    }
    return total;
}

void ChargeDepositor::performDeposit(const std::vector<double>& barycentric,
                                     const std::vector<NodeId>& nodes, double charge_contribution,
                                     std::vector<double>& charge_density)
{
    // 鎸夐噸蹇冨潗鏍囧垎閰嶇數鑽峰埌鍚勪釜鑺傜偣
    for (size_t i = 0; i < 4 && i < nodes.size(); ++i)
    {
        if (nodes[i] < charge_density.size())
        {
            charge_density[nodes[i]] += charge_contribution * barycentric[i];
        }
    }
}

//==============================================================================
// CrossTetrahedron 瀹炵幇
//==============================================================================

CrossTetrahedron::CrossTetrahedron() : boundary_reflection_(false), reflection_potential_(0.0) {}

void CrossTetrahedron::setMesh(const std::vector<Mesh::NodePtr>& nodes,
                               const std::vector<Mesh::ElementPtr>& elements)
{
    nodes_ = nodes;
    elements_ = elements;
}

void CrossTetrahedron::setFieldInterpolator(FieldInterpolatorPtr interpolator)
{
    field_interpolator_ = interpolator;
}

void CrossTetrahedron::setChargeDepositor(ChargeDepositorPtr depositor)
{
    charge_depositor_ = depositor;
}

void CrossTetrahedron::setBoundaryReflection(bool enable, double reflection_potential)
{
    boundary_reflection_ = enable;
    reflection_potential_ = reflection_potential;
}

CrossTetrahedron::TrajectoryResult
CrossTetrahedron::crossTetrahedron(ParticlePtr particle, double dt,
                                   std::vector<double>& charge_density,
                                   bool /* integrate_deposit */, ElementId element_hint)
{
    if (!particle || elements_.empty() || !field_interpolator_)
    {
        TrajectoryResult result;
        result.lost_particle = true;
        result.loss_reason = "Invalid input parameters";
        return result;
    }

    // 璋冪敤 integrateTrajectory锛屽唴閮ㄥ凡閫氳繃 locateParticle 姝ｇ‘瀹氫綅绮掑瓙鎵€鍦ㄥ崟锟?
    return integrateTrajectory(particle, dt, charge_density, false, Utils::Vector3D(),
                               Utils::Vector3D(), element_hint);
}

CrossTetrahedron::TrajectoryResult
CrossTetrahedron::crossTetraUniformE(ParticlePtr particle, double dt,
                                     const Utils::Vector3D& uniform_E,
                                     std::vector<double>& charge_density)
{
    return integrateTrajectory(particle, dt, charge_density, true, uniform_E, Utils::Vector3D(),
                               static_cast<ElementId>(-1));
}

CrossTetrahedron::TrajectoryResult
CrossTetrahedron::crossTetraVaryingE(ParticlePtr particle, double dt,
                                     std::vector<double>& charge_density)
{
    return integrateTrajectory(particle, dt, charge_density, false, Utils::Vector3D(),
                               Utils::Vector3D(), static_cast<ElementId>(-1));
}

CrossTetrahedron::TrajectoryResult CrossTetrahedron::crossTetraUniformEwithB(
    ParticlePtr particle, double dt, const Utils::Vector3D& uniform_E,
    const Utils::Vector3D& uniform_B, std::vector<double>& charge_density)
{
    return integrateTrajectory(particle, dt, charge_density, true, uniform_E, uniform_B,
                               static_cast<ElementId>(-1));
}

CrossTetrahedron::TrajectoryResult CrossTetrahedron::integrateTrajectory(
    ParticlePtr particle, double dt, std::vector<double>& charge_density, bool uniform_field,
    const Utils::Vector3D& E_uniform, const Utils::Vector3D& B_uniform, ElementId element_hint)
{
    TrajectoryResult result;
    result.final_position = particle->getPosition();
    result.final_velocity = particle->getVelocity();
    result.final_element = 0; // 榛樿浣跨敤绗竴涓厓锟?
    result.integration_time = 0.0;
    result.crossed_boundary = false;
    result.lost_particle = false;

    Utils::Point3D current_pos = particle->getPosition();
    Utils::Vector3D current_vel = particle->getVelocity();

    // 淇锛氱敤 locateParticle 鏌ユ壘绮掑瓙鐪熸鎵€鍦ㄧ殑鍗曞厓锛屾秷锟?current_element=0 纭紪锟?bug锟?
    // 鍘熶唬鐮佹案杩滄妸褰撳墠鍗曞厓璁句负 element_0锛堝簳閮ㄥ崟鍏冿級锛屽鑷存敞鍏ヤ簬 z=Lz 鐨勭矑瀛愮殑
    // dichotomyBoundaryIntersection 璧风偣/缁堢偣鍧囦笉锟?element_0 锟?锟?in1==in2==false
    // 锟?姘歌繙杩斿洖 false 锟?lost_particle 涓嶇疆锟?锟?绮掑瓙椋炲嚭鍩熶絾浠嶄负 active锟?
    const bool debug_cross = (std::getenv("PIC_DEBUG_CROSS") != nullptr);

    ElementId current_element = locateParticle(current_pos, element_hint);
    if (current_element == static_cast<ElementId>(-1))
    {
        static std::atomic<int> debug_count{0};
        if (debug_cross)
        {
            int count = debug_count.fetch_add(1);
            if (count < 5)
            {
#pragma omp critical
                {
                    std::cout << "DEBUG locateParticle FAILED at START: pos=(" << current_pos.x()
                              << ", " << current_pos.y() << ", " << current_pos.z()
                              << ") particle charge=" << particle->getCharge()
                              << " mass=" << particle->getMass() << std::endl;
                }
            }
        }

        // 绮掑瓙鍒濆浣嶇疆鍦ㄧ綉鏍间箣澶栵紝鐩存帴鏍囪涓轰涪锟?
        result.lost_particle = true;
        result.loss_reason = "Particle starts outside computational domain";
        return result;
    }

    double remaining_time = dt;
    const double min_time_step = 1e-20;
    const int max_iterations = 100000;
    int iteration = 0;

    double charge = particle->getCharge();
    double mass = particle->getMass();
    double charge_over_mass = charge / mass;

    std::vector<double> barycentric(4); // Hoisted allocation

    while (remaining_time > min_time_step && iteration < max_iterations)
    {
        iteration++;

        // 鑾峰彇褰撳墠浣嶇疆鐨勫満
        Utils::Vector3D E, B;
        if (uniform_field)
        {
            E = E_uniform;
            B = B_uniform;
        }
        else
        {
            E = field_interpolator_->interpolateElectricField(current_pos, current_element);
            B = field_interpolator_->interpolateMagneticField(current_pos, current_element);
        }

        // 璁＄畻鑷€傚簲鏃堕棿姝ラ暱
        double adaptive_dt = adaptiveTimeStep(particle, remaining_time, E, B);

        // 淇濆瓨璧峰浣嶇疆
        Utils::Point3D start_pos = current_pos;

        // 鎵ц绮掑瓙鎺ㄨ繘锛堝畬鏁寸殑Boris绠楁硶锟?
        if (std::abs(charge_over_mass) > 1e-20)
        {
            // 鏈夌數鑽风殑绮掑瓙锛屼娇鐢ㄥ畬鏁寸殑Boris绠楁硶
            current_vel = borisIntegration(current_vel, E, B, charge_over_mass, adaptive_dt);
            current_pos = current_pos + current_vel * adaptive_dt;
        }
        else
        {
            // 涓€х矑瀛愶紝鍖€閫熺洿绾胯繍锟?
            current_pos = current_pos + current_vel * adaptive_dt;
        }

        // 妫€鏌ユ槸鍚︿粛鍦ㄥ綋鍓嶅洓闈綋锟?
        bool in_tetrahedron = field_interpolator_->computeBarycentricCoordinates(
            current_pos, current_element, barycentric);

        if (!in_tetrahedron)
        {
            // 绮掑瓙绂诲紑浜嗗綋鍓嶅洓闈綋锛岄渶瑕佹壘鍒拌竟鐣屼氦锟?
            Utils::Point3D intersection;
            if (detectBoundaryCrossing(start_pos, current_pos, current_element, intersection))
            {
                result.crossed_boundary = true;

                // 澶勭悊杈圭晫鏉′欢
                if (boundary_reflection_ && handleBoundaryCondition(particle, intersection))
                {
                    // 鍙嶅皠澶勭悊
                    current_pos = intersection;
                    stats_.reflection_events++;
                }
                else
                {
                    // 灏濊瘯鎵惧埌鏂扮殑鍥涢潰锟?
                    // Reuse previous element as hint to trigger neighbor-first search.
                    ElementId new_element = locateParticle(current_pos, current_element);
                    if (new_element != static_cast<ElementId>(-1))
                    {
                        current_element = new_element;
                    }
                    else
                    {
                        static std::atomic<int> lost_debug_count{0};
                        if (debug_cross)
                        {
                            int count = lost_debug_count.fetch_add(1);
                            if (count < 5)
                            {
#pragma omp critical
                                {
                                    std::cout
                                        << "DEBUG locateParticle FAILED after MOVE: start_pos=("
                                        << start_pos.x() << ", " << start_pos.y() << ", "
                                        << start_pos.z() << ") end_pos=(" << current_pos.x() << ", "
                                        << current_pos.y() << ", " << current_pos.z()
                                        << ") old_elem=" << current_element << std::endl;
                                }
                            }
                        }

                        // 绮掑瓙涓㈠け
                        result.lost_particle = true;
                        result.loss_reason = "Particle left computational domain";
                        break;
                    }
                }
            }
        }

        // 娌夌Н鐢佃嵎
        if (charge_depositor_ && !charge_density.empty())
        {
            charge_depositor_->depositChargeAlongTrajectory(start_pos, current_pos, current_element,
                                                            charge, particle->getWeight(),
                                                            adaptive_dt, charge_density);
        }

        remaining_time -= adaptive_dt;
        result.integration_time += adaptive_dt;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
        stats_.total_crossings++;
    }

    if (iteration >= max_iterations)
    {
        static std::atomic<int> iter_debug_count{0};
        if (debug_cross)
        {
            int count = iter_debug_count.fetch_add(1);
            if (count < 5)
            {
#pragma omp critical
                {
                    std::cout << "DEBUG MAX ITERATIONS EXCEEDED: pos=(" << current_pos.x() << ", "
                              << current_pos.y() << ", " << current_pos.z()
                              << ") iterations=" << iteration << std::endl;
                }
            }
        }

        result.lost_particle = true;
        result.loss_reason = "Maximum iterations exceeded";

#ifdef USE_OPENMP
#pragma omp atomic
#endif
        stats_.lost_particles++;
    }

    // 鏇存柊缁撴灉
    result.final_position = current_pos;
    result.final_velocity = current_vel;
    result.final_element = current_element;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
    stats_.total_integration_time += result.integration_time;

    return result;
}

bool CrossTetrahedron::detectBoundaryCrossing(const Utils::Point3D& start,
                                              const Utils::Point3D& end, ElementId element_id,
                                              Utils::Point3D& intersection)
{
    // 浣跨敤浜屽垎娉曟壘鍒拌竟鐣屼氦锟?
    return dichotomyBoundaryIntersection(start, end, element_id, intersection);
}

double CrossTetrahedron::adaptiveTimeStep(ParticlePtr particle, double dt_max,
                                          const Utils::Vector3D& E, const Utils::Vector3D& B)
{
    // 绠€鍖栫殑鑷€傚簲鏃堕棿姝ラ暱绠楁硶
    double dt = dt_max;

    // 鍩轰簬鐢靛満寮哄害闄愬埗鏃堕棿姝ラ暱
    double E_magnitude = E.magnitude();
    if (E_magnitude > 1e-10)
    {
        double charge_over_mass = particle->getCharge() / particle->getMass();
        double acceleration = E_magnitude * std::abs(charge_over_mass);
        double velocity = particle->getVelocity().magnitude();

        // 闄愬埗閫熷害鍙樺寲涓嶈秴锟?0%
        if (acceleration > 1e-10 && velocity > 1.0)
        { // 娣诲姞velocity > 1m/s鐨勪繚锟?
            double dt_velocity = 0.1 * velocity / acceleration;
            dt = std::min(dt, dt_velocity);
        }
    }

    // 鍩轰簬纾佸満寮哄害闄愬埗鏃堕棿姝ラ暱锛堝洖鏃嬮鐜囷級
    double B_magnitude = B.magnitude();
    if (B_magnitude > 1e-10)
    {
        double charge_over_mass = particle->getCharge() / particle->getMass();
        double cyclotron_freq = B_magnitude * std::abs(charge_over_mass);
        if (cyclotron_freq > 1e-10)
        {
            double dt_cyclotron = 0.1 * 2.0 * M_PI / cyclotron_freq;
            dt = std::min(dt, dt_cyclotron);
        }
    }

    // 纭繚鏃堕棿姝ラ暱鍦ㄥ悎鐞嗚寖鍥村唴锟?
    // 榛樿 floor_ratio=0.001 淇濇寔鍘熸湁鐗╃悊琛屼负锟?
    // 鍙€氳繃鐜鍙橀噺 PIC_ADAPTIVE_FLOOR_RATIO 涓存椂璋冨ぇ浠ュ仛鎬ц兘娴嬭瘯锟?
    double floor_ratio = 0.001;
    if (const char* env_ratio = std::getenv("PIC_ADAPTIVE_FLOOR_RATIO"))
    {
        const double parsed = std::atof(env_ratio);
        if (std::isfinite(parsed) && parsed > 0.0)
        {
            floor_ratio = std::min(std::max(parsed, 1.0e-5), 0.5);
        }
    }
    dt = std::max(dt, dt_max * floor_ratio);
    dt = std::min(dt, dt_max);

    return dt;
}

Utils::Vector3D CrossTetrahedron::borisIntegration(const Utils::Vector3D& velocity,
                                                   const Utils::Vector3D& E,
                                                   const Utils::Vector3D& B,
                                                   double charge_over_mass, double dt)
{
    /**
     * Boris绠楁硶鏄疨IC妯℃嫙涓殑鏍囧噯绮掑瓙鎺ㄨ繘绠楁硶
     * 瀹冭兘澶熺簿纭鐞嗗甫鐢电矑瀛愬湪鐢电鍦轰腑鐨勮繍锟?
     *
     * 绠楁硶姝ラ锟?
     * 1. v- = v + (q/m) * E * dt/2  (鐢靛満鍔犻€熶竴鍗婃椂闂存)
     * 2. 纾佸満鏃嬭浆
     * 3. v+ = v- + (q/m) * E * dt/2  (鐢靛満鍔犻€熷彟涓€鍗婃椂闂存)
     */

    const double half_dt = 0.5 * dt;

    // 姝ラ1: 鐢靛満鍔犻€熶竴鍗婃椂闂存
    Utils::Vector3D v_minus = velocity + E * (charge_over_mass * half_dt);

    // 姝ラ2: 纾佸満鏃嬭浆
    Utils::Vector3D v_prime = v_minus;

    double B_magnitude = B.magnitude();
    if (B_magnitude > 1e-20)
    {
        // 璁＄畻鏃嬭浆鍙傛暟
        Utils::Vector3D t = B * (charge_over_mass * half_dt);
        double t_magnitude_squared = t.dot(t);

        // 閬垮厤鏁板€间笉绋冲畾
        if (t_magnitude_squared < 100.0)
        { // |蠅t| < 10
            Utils::Vector3D s = t * (2.0 / (1.0 + t_magnitude_squared));

            // Boris鏃嬭浆
            Utils::Vector3D v_cross_t = v_minus.cross(t);
            Utils::Vector3D v_star = v_minus + v_cross_t;
            Utils::Vector3D v_star_cross_s = v_star.cross(s);
            v_prime = v_minus + v_star_cross_s;
        }
        // 濡傛灉纾佸満澶己锛岃烦杩囨棆杞楠や互閬垮厤鏁板€间笉绋冲畾
    }

    // 姝ラ3: 鐢靛満鍔犻€熷彟涓€鍗婃椂闂存
    Utils::Vector3D v_plus = v_prime + E * (charge_over_mass * half_dt);

    return v_plus;
}

ElementId CrossTetrahedron::locateParticle(const Utils::Point3D& position, ElementId hint)
{
    if (elements_.empty())
        return static_cast<ElementId>(-1);

    std::vector<double> barycentric(4); // Hoisted allocation

    // 1. 棣栧厛妫€鏌ユ彁绀哄崟锟?
    if (hint != static_cast<ElementId>(-1) && hint < static_cast<ElementId>(elements_.size()))
    {
        if (field_interpolator_->computeBarycentricCoordinates(position, hint, barycentric))
        {
            return hint;
        }

        // 2. 妫€鏌ユ彁绀哄崟鍏冪殑閭诲眳 (Optimization: Neighbor Search)
        const auto& element = elements_[hint];
        const auto& neighbors = element->getNeighbors();
        for (const auto& neighbor : neighbors)
        {
            if (field_interpolator_->computeBarycentricCoordinates(position, neighbor->getId(),
                                                                   barycentric))
            {
                return neighbor->getId();
            }
        }
    }

    // 3. 閬嶅巻鎵€鏈夊崟鍏冿紙浣滀负鏈€鍚庣殑鎵嬫锟?
    for (ElementId i = 0; i < static_cast<ElementId>(elements_.size()); ++i)
    {
        if (field_interpolator_->computeBarycentricCoordinates(position, i, barycentric))
        {
            return i;
        }
    }

    return static_cast<ElementId>(-1); // 鏈壘锟?
}

bool CrossTetrahedron::dichotomyBoundaryIntersection(const Utils::Point3D& start,
                                                     const Utils::Point3D& end,
                                                     ElementId element_id,
                                                     Utils::Point3D& intersection, double tolerance)
{
    Utils::Point3D p1 = start;
    Utils::Point3D p2 = end;

    // 妫€鏌ヨ捣鐐规槸鍚﹀湪鍥涢潰浣撳唴
    std::vector<double> barycentric1(4), barycentric2(4);
    bool in1 = field_interpolator_->computeBarycentricCoordinates(p1, element_id, barycentric1);
    bool in2 = field_interpolator_->computeBarycentricCoordinates(p2, element_id, barycentric2);

    if (in1 == in2)
    {
        return false; // 娌℃湁杈圭晫绌胯秺
    }

    std::vector<double> barycentric_mid(4); // Hoisted allocation

    // 浜屽垎娉曟煡鎵捐竟鐣屼氦锟?
    for (int i = 0; i < 50; ++i)
    { // 鏈€锟?0娆¤凯锟?
        Utils::Point3D mid = (p1 + p2) * 0.5;

        if ((p2 - p1).magnitude() < tolerance)
        {
            intersection = mid;
            return true;
        }

        bool in_mid =
            field_interpolator_->computeBarycentricCoordinates(mid, element_id, barycentric_mid);

        if (in_mid == in1)
        {
            p1 = mid;
        }
        else
        {
            p2 = mid;
        }
    }

    intersection = (p1 + p2) * 0.5;
    return true;
}

bool CrossTetrahedron::handleBoundaryCondition(ParticlePtr particle,
                                               const Utils::Point3D& intersection_point)
{
    if (!boundary_reflection_)
    {
        return false;
    }

    // 绠€鍖栫殑鍙嶅皠澶勭悊
    // 瀹為檯搴旂敤涓渶瑕佹牴鎹竟鐣屾硶鍚戦噺璁＄畻鍙嶅皠閫熷害
    Utils::Vector3D velocity = particle->getVelocity();

    // 鍋囪杈圭晫娉曞悜閲忎负z鏂瑰悜锛堢畝鍖栵級
    Utils::Vector3D normal(0.0, 0.0, 1.0);

    // 璁＄畻鍙嶅皠閫熷害
    Utils::Vector3D reflected_velocity = velocity - normal * (2.0 * velocity.dot(normal));

    particle->setVelocity(reflected_velocity);
    particle->setPosition(intersection_point);

    return true;
}

} // namespace Particle
} // namespace SCDAT

