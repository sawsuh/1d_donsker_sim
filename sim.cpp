#include <cmath>
#include <iostream>
#include <unordered_map>

enum PlusMinus { plus, minus };

float a(float x) { return 1; }

float rho(float x) { return 1; }

float b(float x) { return 0; }

struct cellData {
  float time_left;
  float time_right;
  float prob_left;
  float prob_right;
};

class CellDataCalculator {
public:
  float left;
  float right;
  float x;
  CellDataCalculator(float l, float r, float y) {
    left = l;
    right = r;
    x = y;
    integration_inc = 0.001;
  }
  cellData compute_cell_data() {
    float prob_left = v0minus(x);
    float prob_right = v0plus(x);
    float time_left_ind = v1minus(x);
    float time_right_ind = v1plus(x);
    return cellData{time_left_ind / prob_left, time_right_ind / prob_right,
                    prob_left, prob_right};
  }

private:
  std::unordered_map<float, float> psi_values;
  std::unordered_map<float, float> v0plus_helper_values;
  // std::unordered_map<float, float> v0plus_values;
  // std::unordered_map<float, float> v0minus_values;
  float integration_inc;
  float psi(float x) {
    if (psi_values.find(x) != psi_values.end()) {
      return psi_values[x];
    }
    float integral = 0;
    float y = left;
    while (y < right) {
      integral += b(y) * integration_inc / (a(y) * rho(y));
      y += integration_inc;
    }
    integral *= 2;
    psi_values[x] = integral;
    return integral;
  }
  float v0plus_helper_integral(float x) {
    if (v0plus_helper_values.find(x) != v0plus_helper_values.end()) {
      return v0plus_helper_values[x];
    }
    float integral = 0;
    float y = left;
    while (y < right) {
      integral += integration_inc * exp(-psi(y)) / a(y);
      y += integration_inc;
    }
    v0plus_helper_values[x] = integral;
    return integral;
  }
  float v0plus(float x) {
    // if (v0plus_values.find(x) != v0plus_values.end()) {
    //   return v0plus_values[x];
    // }
    // v0plus_values[x] =
    // v0plus_helper_integral(x) / v0plus_helper_integral(right);
    // return v0plus_values[x];
    return v0plus_helper_integral(x) / v0plus_helper_integral(right);
  }
  float v0minus(float x) {
    // if (v0minus_values.find(x) != v0minus_values.end()) {
    //   return v0minus_values[x];
    // }
    // v0minus_values[x] = 1 - v0plus(x);
    // return v0minus_values[x];
    return 1 - v0plus(x);
  }
  float v1plus(float x) { return 0; }
  float v1minus(float x) { return 0; }
};

int main() {
  std::cout << "hi" << std::endl;
  return 0;
}
