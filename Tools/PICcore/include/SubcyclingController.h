#ifndef SCDAT_PICCORE_SUBCYCLINGCONTROLLER_H
#define SCDAT_PICCORE_SUBCYCLINGCONTROLLER_H

/**
 * @file SubcyclingController.h
 * @brief Subcycling time-step controller for electron and ion updates.
 */

#include <cmath>
#include <stdexcept>
#include <string>

namespace SCDAT
{
namespace PICcore
{

class SubcyclingController
{
  public:
    struct Parameters
    {
        double electron_dt = 1e-12;
        int subcycle_ratio = 20;
        bool enabled = true;

        void validate() const
        {
            if (electron_dt <= 0.0)
            {
                throw std::invalid_argument("electron_dt must be positive");
            }
            if (subcycle_ratio < 1)
            {
                throw std::invalid_argument("subcycle_ratio must be >= 1");
            }
            if (subcycle_ratio > 100)
            {
                throw std::invalid_argument(
                    "subcycle_ratio too large (> 100), may cause instability");
            }
        }

        double getIonTimeStep() const
        {
            return subcycle_ratio * electron_dt;
        }
    };

    struct Statistics
    {
        int total_electron_steps = 0;
        int total_ion_steps = 0;
        double total_time = 0.0;

        double getSpeedup() const
        {
            if (total_electron_steps == 0 || total_ion_steps == 0)
            {
                return 1.0;
            }
            return static_cast<double>(total_electron_steps) / static_cast<double>(total_ion_steps);
        }

        void reset()
        {
            total_electron_steps = 0;
            total_ion_steps = 0;
            total_time = 0.0;
        }
    };

    explicit SubcyclingController(int subcycle_ratio = 20, double electron_dt = 1e-12) : counter_(0)
    {
        params_.subcycle_ratio = subcycle_ratio;
        params_.electron_dt = electron_dt;
        params_.enabled = true;
        params_.validate();
    }

    explicit SubcyclingController(const Parameters& params) : params_(params), counter_(0)
    {
        params_.validate();
    }

    bool shouldPushIons() const
    {
        if (!params_.enabled)
        {
            return true;
        }
        return (counter_ % params_.subcycle_ratio) == 0;
    }

    void incrementCounter()
    {
        counter_++;
        stats_.total_electron_steps++;
        stats_.total_time += params_.electron_dt;

        if (shouldPushIons())
        {
            stats_.total_ion_steps++;
        }
    }

    void reset(bool reset_stats = false)
    {
        counter_ = 0;
        if (reset_stats)
        {
            stats_.reset();
        }
    }

    double getElectronTimeStep() const
    {
        return params_.electron_dt;
    }
    double getIonTimeStep() const
    {
        return params_.getIonTimeStep();
    }
    int getCounter() const
    {
        return counter_;
    }
    int getSubcycleRatio() const
    {
        return params_.subcycle_ratio;
    }

    void setSubcycleRatio(int ratio)
    {
        if (ratio < 1 || ratio > 100)
        {
            throw std::invalid_argument("Invalid subcycle_ratio: must be in [1, 100]");
        }
        params_.subcycle_ratio = ratio;
        reset();
    }

    void setElectronTimeStep(double dt)
    {
        if (dt <= 0.0)
        {
            throw std::invalid_argument("electron_dt must be positive");
        }
        params_.electron_dt = dt;
    }

    void setEnabled(bool enabled)
    {
        params_.enabled = enabled;
        if (!enabled)
        {
            reset();
        }
    }

    bool isEnabled() const
    {
        return params_.enabled;
    }
    const Parameters& getParameters() const
    {
        return params_;
    }
    const Statistics& getStatistics() const
    {
        return stats_;
    }
    double getSpeedup() const
    {
        return stats_.getSpeedup();
    }

    std::string getStatusString() const
    {
        std::string status;
        status += "SubcyclingController Status:\n";
        status += "  Enabled: " + std::string(params_.enabled ? "YES" : "NO") + "\n";
        status += "  Electron dt: " + std::to_string(params_.electron_dt) + " s\n";
        status += "  Ion dt: " + std::to_string(getIonTimeStep()) + " s\n";
        status += "  Subcycle ratio: " + std::to_string(params_.subcycle_ratio) + "\n";
        status += "  Current counter: " + std::to_string(counter_) + "\n";
        status += "  Total electron steps: " + std::to_string(stats_.total_electron_steps) + "\n";
        status += "  Total ion steps: " + std::to_string(stats_.total_ion_steps) + "\n";
        status += "  Total time: " + std::to_string(stats_.total_time) + " s\n";
        status += "  Speedup: " + std::to_string(getSpeedup()) + "x\n";
        return status;
    }

  private:
    Parameters params_;
    Statistics stats_;
    int counter_;
};

} // namespace PICcore
} // namespace SCDAT

#endif // SCDAT_PICCORE_SUBCYCLINGCONTROLLER_H
