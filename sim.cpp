#include <iostream>
#include <map>

enum PlusMinus { plus, minus };

float a(float x) { return 1; }

float rho(float x) { return 1; }

float b(float x) { return 0; }

float psi(float x0, float x, std::map<float, float> cache) { return 2; }

struct cellData {
  float time_left;
  float time_right;
  float prob_left;
  float prob_right;
};

cellData compute_cell_data(float x, float left, float right,
                           float integration_inc = 0.001) {
  std::map<float, float> psi_values;
  std::map<float, float> v0plus_values;
  std::map<float, float> v0minus_values;
  // return 2 * b(x) / (a(x) * rho(x));
  float prob_left = psi(left, x, psi_values);
  return cellData{0, 0, prob_left, 0};
}

int main() {
  std::cout << "hi" << std::endl;
  return 0;
}
