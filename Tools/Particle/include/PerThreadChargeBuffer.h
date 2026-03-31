#ifndef SCDAT_PER_THREAD_CHARGE_BUFFER_H
#define SCDAT_PER_THREAD_CHARGE_BUFFER_H

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

namespace SCDAT
{
namespace Particle
{

/**
 * @brief Thread-local charge accumulation buffers for parallel deposition.
 */
class PerThreadChargeBuffer
{
  public:
    PerThreadChargeBuffer() : node_count_(0), thread_count_(1) {}

    explicit PerThreadChargeBuffer(std::size_t node_count, int thread_count = 0)
    {
        resize(node_count, thread_count);
    }

    void resize(std::size_t node_count, int thread_count = 0)
    {
        node_count_ = node_count;

        if (thread_count <= 0)
        {
            const unsigned int hardware_threads = std::thread::hardware_concurrency();
            thread_count_ = std::max(1, static_cast<int>(hardware_threads == 0 ? 1 : hardware_threads));
        }
        else
        {
            thread_count_ = std::max(1, thread_count);
        }

        buffers_.assign(static_cast<std::size_t>(thread_count_),
                        std::vector<double>(node_count_, 0.0));
    }

    void clear()
    {
        for (auto& buffer : buffers_)
        {
            std::fill(buffer.begin(), buffer.end(), 0.0);
        }
    }

    std::vector<double>& localBuffer(int thread_index)
    {
        return buffers_.at(normalizeThreadIndex(thread_index));
    }

    const std::vector<double>& localBuffer(int thread_index) const
    {
        return buffers_.at(normalizeThreadIndex(thread_index));
    }

    void addCharge(int thread_index, std::size_t node_index, double charge)
    {
        buffers_.at(normalizeThreadIndex(thread_index)).at(node_index) += charge;
    }

    std::vector<double> reduce() const
    {
        std::vector<double> reduced(node_count_, 0.0);
        reduceInto(reduced);
        return reduced;
    }

    void reduceInto(std::vector<double>& target) const
    {
        target.assign(node_count_, 0.0);
        for (const auto& buffer : buffers_)
        {
            for (std::size_t i = 0; i < node_count_; ++i)
            {
                target[i] += buffer[i];
            }
        }
    }

    std::size_t nodeCount() const { return node_count_; }
    int threadCount() const { return thread_count_; }

  private:
    std::size_t normalizeThreadIndex(int thread_index) const
    {
        if (buffers_.empty())
        {
            return 0;
        }

        if (thread_index < 0)
        {
            return 0;
        }

        return std::min<std::size_t>(static_cast<std::size_t>(thread_index), buffers_.size() - 1);
    }

  private:
    std::size_t node_count_;
    int thread_count_;
    std::vector<std::vector<double>> buffers_;
};

} // namespace Particle
} // namespace SCDAT

#endif // SCDAT_PER_THREAD_CHARGE_BUFFER_H
