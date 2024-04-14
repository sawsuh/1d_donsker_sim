#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <unistd.h>
#include <unordered_map>

enum PlusMinus { plus, minus };

double a(double x) { return 1; }

double rho(double x) { return 1; }

double b(double x) { return 0; }

struct cellData {
  double time_left;
  double time_right;
  double prob_left;
  double prob_right;
};

class CellDataCalculator {
public:
  double left;
  double right;
  double x;
  CellDataCalculator(double l, double r, double y) {
    left = l;
    right = r;
    x = y;
    integration_inc = 0.00001;
  }
  cellData compute_cell_data() {
    gen_psi_table();
    gen_v0plus_helper_table();
    double prob_left = v0minus(x);
    double prob_right = v0plus(x);
    double time_left_ind = v1minus(x);
    double time_right_ind = v1plus(x);
    return cellData{time_left_ind / prob_left, time_right_ind / prob_right,
                    prob_left, prob_right};
  }

private:
  std::unordered_map<double, double> psi_values;
  std::unordered_map<double, double> v0plus_helper_values;
  double integration_inc;
  void gen_psi_table() {
    double integral = 0;
    double y = left;
    while (y < right) {
      double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (b(y) / (a(y) * rho(y)) + b(y_next) / (a(y_next) * rho(y_next))) / 2;
      psi_values[y] = integral;
      y = y_next;
    }
  }
  double psi(double x) {
    if (psi_values.find(x) != psi_values.end()) {
      return psi_values[x];
    }
    std::cout << "psi cache miss" << std::endl;
    double integral = 0;
    double y = left;
    while (y < x) {
      double y_next = y + integration_inc;
      if (psi_values.find(y) != psi_values.end()) {
        integral = psi_values[y];
        y = y_next;
        continue;
      }
      integral +=
          integration_inc *
          (b(y) / (a(y) * rho(y)) + b(y_next) / (a(y_next) * rho(y_next))) / 2;
      psi_values[y] = integral;
      y = y_next;
    }
    integral *= 2;
    psi_values[x] = integral;
    return integral;
  }
  void gen_v0plus_helper_table() {
    double integral = 0;
    double y = left;
    while (y < right) {
      double y_next = y + integration_inc;
      integral += integration_inc *
                  (exp(-psi(y)) / a(y) + exp(-psi(y_next)) / a(y_next)) / 2;
      v0plus_helper_values[y] = integral;
      y = y_next;
    }
  }
  double v0plus_helper_integral(double x) {
    if (v0plus_helper_values.find(x) != v0plus_helper_values.end()) {
      return v0plus_helper_values[x];
    }
    double integral = 0;
    double y = left;
    std::cout << "v0plus helper cache miss" << std::endl;
    while (y < x) {
      double y_next = y + integration_inc;
      if (v0plus_helper_values.find(y) != v0plus_helper_values.end()) {
        integral = v0plus_helper_values[y];
        y = y_next;
        continue;
      }
      integral += integration_inc *
                  (exp(-psi(y)) / a(y) + exp(-psi(y_next)) / a(y_next)) / 2;
      v0plus_helper_values[y] = integral;
      y = y_next;
    }
    v0plus_helper_values[x] = integral;
    return integral;
  }
  double v0plus(double x) {
    return v0plus_helper_integral(x) / v0plus_helper_integral(right);
  }
  double v0minus(double x) { return 1 - v0plus(x); }
  double G(double x, double y) {
    if (x <= y) {
      return 2 * (v0plus_helper_integral(x) - v0plus_helper_integral(left)) *
             (v0plus_helper_integral(right) - v0plus_helper_integral(y)) /
             (v0plus_helper_integral(right) - v0plus_helper_integral(left));
    } else {
      return 2 * (v0plus_helper_integral(y) - v0plus_helper_integral(left)) *
             (v0plus_helper_integral(right) - v0plus_helper_integral(x)) /
             (v0plus_helper_integral(right) - v0plus_helper_integral(left));
    }
  }
  double v1plus(double x) {
    double integral = 0;
    double y = left;
    while (y < right) {
      double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (G(x, y) * v0plus(y) * exp(psi(y)) / rho(y) +
           G(x, y_next) * v0plus(y_next) * exp(psi(y_next)) / rho(y_next)) /
          2;
      y = y_next;
    }
    return integral;
  }
  double v1minus(double x) {
    double integral = 0;
    double y = left;
    while (y < right) {
      double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (G(x, y) * v0minus(y) * exp(psi(y)) / rho(y) +
           G(x, y_next) * v0minus(y_next) * exp(psi(y_next)) / rho(y_next)) /
          2;
      y = y_next;
    }
    return integral;
  }
};

struct cell {
  double left, right;
};
struct increment {
  double next_point, delta_t;
};
class Simulator {
public:
  double start;
  Simulator(double x) { start = x; }
  void simulate(double t, int rounds = 1000) {
    std::vector<double> out;
    for (auto _ = rounds; _--;) {
      double cur = start;
      double t_cur = 0;

      while (t_cur < t) {
        increment inc = next_point(cur);
        cur = inc.next_point;
        t_cur += inc.delta_t;
      }
      std::cout << cur << std::endl;
      out.push_back(cur);
    }
    std::ofstream output_file("res.csv");
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    std::copy(std::begin(out), std::end(out), output_iterator);
  }

private:
  std::unordered_map<double, cellData> cell_cache;
  std::random_device rd;
  std::mt19937 rng{rd()};
  // INPUT grid spec
  cell get_adjacent(double point) { return cell{point - 0.05, point + 0.05}; }
  cellData get_data(double point) {
    if (cell_cache.find(point) != cell_cache.end()) {
      return cell_cache[point];
    }
    std::cout << "cell cache miss" << std::endl;
    cell lr = get_adjacent(point);
    CellDataCalculator calc(lr.left, lr.right, point);
    cellData out = calc.compute_cell_data();
    cell_cache[point] = out;
    return out;
  }
  increment next_point(double point) {
    cell lr = get_adjacent(point);
    cellData point_data = get_data(point);
    std::bernoulli_distribution d(point_data.prob_right);
    if (d(rng)) { // exit right
      return increment{lr.right, point_data.time_right};
    } else { // exit left
      return increment{lr.left, point_data.time_left};
    }
  }
};

int main() {
  // INPUT
  double start = 0;
  double time = 1;
  Simulator sim(start);
  sim.simulate(time);
  return 0;
}
