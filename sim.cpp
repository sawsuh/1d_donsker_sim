#include <cmath>
#include <iostream>
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
    integration_inc = 0.001;
  }
  cellData compute_cell_data() {
    double prob_left = v0minus(x);
    double prob_right = v0plus(x);
    double time_left_ind = v1minus(x);
    double time_right_ind = v1plus(x);
    return cellData{time_left_ind / prob_left, time_right_ind / prob_right,
                    prob_left, prob_right};
  }
  double test() {
    double integral = 0;
    double y = left;
    while (y < right) {
      if (x <= y) {
        integral +=
            integration_inc * 2 * (x - left) * (right - y) / (right - left);
      } else {
        integral +=
            integration_inc * 2 * (y - left) * (right - x) / (right - left);
      }
      y += integration_inc;
    }
    v0plus_helper_values[x] = integral;
    return integral;
  }

private:
  std::unordered_map<double, double> psi_values;
  std::unordered_map<double, double> v0plus_helper_values;
  std::unordered_map<double, double> v1plus_values;
  std::unordered_map<double, double> v1minus_values;
  // std::unordered_map<double, double> v0plus_values;
  // std::unordered_map<double, double> v0minus_values;
  double integration_inc;
  double psi(double x) {
    if (psi_values.find(x) != psi_values.end()) {
      // std::cout << "psi(" << x << ") = " << psi_values[x] << std::endl;
      return psi_values[x];
    }
    double integral = 0;
    double y = left;
    while (y < right) {
      double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (b(y) / (a(y) * rho(y)) + b(y_next) / (a(y_next) * rho(y_next))) / 2;
      y = y_next;
    }
    integral *= 2;
    psi_values[x] = integral;
    return integral;
  }
  double v0plus_helper_integral(double x) {
    if (v0plus_helper_values.find(x) != v0plus_helper_values.end()) {
      // std::cout << "v0pint(" << x << ") = " << v0plus_helper_values[x]
      //           << std::endl;
      return v0plus_helper_values[x];
    }
    double integral = 0;
    double y = left;
    while (y < x) {
      double y_next = y + integration_inc;
      integral += integration_inc *
                  (exp(-psi(y)) / a(y) + exp(-psi(y_next)) / a(y_next)) / 2;
      y = y_next;
    }
    v0plus_helper_values[x] = integral;
    return integral;
  }
  double v0plus(double x) {
    // if (v0plus_values.find(x) != v0plus_values.end()) {
    //   return v0plus_values[x];
    // }
    // v0plus_values[x] =
    // v0plus_helper_integral(x) / v0plus_helper_integral(right);
    // return v0plus_values[x];
    return v0plus_helper_integral(x) / v0plus_helper_integral(right);
  }
  double v0minus(double x) {
    // if (v0minus_values.find(x) != v0minus_values.end()) {
    //   return v0minus_values[x];
    // }
    // v0minus_values[x] = 1 - v0plus(x);
    // return v0minus_values[x];
    return 1 - v0plus(x);
  }
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
    if (v1plus_values.find(x) != v1plus_values.end()) {
      // std::cout << "psi(" << x << ") = " << psi_values[x] << std::endl;
      return v1plus_values[x];
    }
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
    v1plus_values[x] = integral;
    return integral;
  }
  double v1minus(double x) {
    if (v1minus_values.find(x) != v1minus_values.end()) {
      // std::cout << "psi(" << x << ") = " << psi_values[x] << std::endl;
      return v1minus_values[x];
    }
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
    v1minus_values[x] = integral;
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
      // out.push_back(cur);
      std::cout << cur << std::endl;
    }
    // return out;
  }

private:
  std::unordered_map<double, cellData> cell_cache;
  std::random_device rd;
  std::mt19937 rng{rd()};
  cell get_adjacent(double point) { return cell{point - 0.1, point + 0.1}; }
  cellData get_data(double point) {
    if (cell_cache.find(point) != cell_cache.end()) {
      return cell_cache[point];
    }
    cell lr = get_adjacent(point);
    CellDataCalculator calc(lr.left, lr.right, point);
    cellData out = calc.compute_cell_data();
    cell_cache[point] = out;
    return out;
  }
  increment next_point(double point) {
    cell lr = get_adjacent(point);
    cellData point_data = get_data(point);
    /* std::cout << "left: " << point_data.prob_left
              << ", right: " << point_data.prob_right << std::endl; */
    std::bernoulli_distribution d(point_data.prob_right);
    // sleep(1);
    if (d(rng)) { // exit right
      /* std::cout << "right from " << point << " to " << lr.right << " in "
                << point_data.time_right << std::endl; */
      return increment{lr.right, point_data.time_right};
    } else { // exit left
      /* std::cout << "left from " << point << " to " << lr.left << "in "
                << point_data.time_left << std::endl; */
      return increment{lr.left, point_data.time_left};
    }
  }
};

int main() {
  /* double l = -0.2;
  double r = 0.3;
  double x = 0;
  CellDataCalculator calc(l, r, x);
  cellData dat = calc.compute_cell_data();
  std::cout << dat.prob_left << std::endl;
  std::cout << dat.prob_right << std::endl;
  std::cout << dat.time_left << std::endl;
  std::cout << dat.time_right << std::endl;
  std::cout << "we computed: "
            << dat.time_right * dat.prob_right + dat.time_left * dat.prob_left
            << std::endl;
  std::cout << "we also computed: " << calc.test() << std::endl;
  std::cout << "true value: " << -l * r << std::endl; */

  Simulator sim(0);
  sim.simulate(1);
  return 0;
}
