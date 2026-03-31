#pragma once

#include <stdexcept>
#include <vector>

namespace SCDAT
{
namespace Core
{

class BandMatrix
{
  public:
    BandMatrix() = default;

    explicit BandMatrix(int size) : values_(static_cast<size_t>(size), std::vector<double>(static_cast<size_t>(size), 0.0))
    {
    }

    BandMatrix(int size, int /*bandwidth*/) : BandMatrix(size) {}

    int getSize() const
    {
        return static_cast<int>(values_.size());
    }

    void resize(int size)
    {
        values_.assign(static_cast<size_t>(size),
                       std::vector<double>(static_cast<size_t>(size), 0.0));
    }

    double& operator()(int row, int col)
    {
        return values_.at(static_cast<size_t>(row)).at(static_cast<size_t>(col));
    }

    const double& operator()(int row, int col) const
    {
        return values_.at(static_cast<size_t>(row)).at(static_cast<size_t>(col));
    }

    void set(int row, int col, double value)
    {
        (*this)(row, col) = value;
    }

    double get(int row, int col) const
    {
        return (*this)(row, col);
    }

    void addValue(int row, int col, double value)
    {
        (*this)(row, col) += value;
    }

    void multiply(const std::vector<double>& input, std::vector<double>& output) const
    {
        const size_t n = values_.size();
        output.assign(n, 0.0);
        for (size_t i = 0; i < n; ++i)
        {
            double sum = 0.0;
            for (size_t j = 0; j < values_[i].size() && j < input.size(); ++j)
            {
                sum += values_[i][j] * input[j];
            }
            output[i] = sum;
        }
    }

  private:
    std::vector<std::vector<double>> values_;
};

} // namespace Core
} // namespace SCDAT
