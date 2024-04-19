#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <stdexcept>
#include <thread>
#include <unistd.h>
#include <unordered_map>

enum PlusMinus { plus, minus };

// INPUT
const int ROUNDS = 10000;
const int PRINT_INTERVAL = 100;
const double INTEGRATION_INC = 0.00001;
const double START = 0;
const double TIME = 1;
double a(double x) { return 1; }
double rho(double x) {
  if (x < 0) {
    return 1;
  } else {
    return 5;
  }
}
double b(double x) { return 0; }
struct cell {
  double left, right;
};
cell get_adjacent(double point) { return cell{point - 0.01, point + 0.01}; }

struct cellData {
  double time_left;
  double time_right;
  double prob_left;
  double prob_right;
};

template <typename DT> void push_safe(std::vector<DT> &vec, int idx, DT x) {
  if (vec.size() != idx) {
    throw std::out_of_range("desired location not end of vector");
  }
  vec.push_back(x);
}

class CellDataCalculator {
public:
  double left;
  double right;
  double x;
  CellDataCalculator(double l, double r, double y) {
    left = l;
    right = r;
    x = y;
    integration_inc = INTEGRATION_INC;
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
  std::vector<double> psi_values;
  std::vector<double> v0plus_helper_values;
  double integration_inc;
  int get_index(double x) { return round((x - left) / integration_inc); }
  void gen_psi_table() {
    psi_values.reserve(get_index(right));
    psi_values.push_back(0);
    double integral = 0;
    double y = left;
    while (y < right) {
      double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (b(y) / (a(y) * rho(y)) + b(y_next) / (a(y_next) * rho(y_next))) / 2;
      push_safe(psi_values, get_index(y_next), integral);
      y = y_next;
    }
  }
  double psi(double x) {

    int position = get_index(x);
    if (position < psi_values.size()) {
      return psi_values[position] * 2;
    }
    throw std::out_of_range("psi cache miss");
  }
  void gen_v0plus_helper_table() {
    v0plus_helper_values.reserve(get_index(right));
    v0plus_helper_values.push_back(0);
    double integral = 0;
    double y = left;
    while (y < right) {
      double y_next = y + integration_inc;
      integral += integration_inc *
                  (exp(-psi(y)) / a(y) + exp(-psi(y_next)) / a(y_next)) / 2;
      push_safe(v0plus_helper_values, get_index(y_next), integral);
      y = y_next;
    }
  }
  double v0plus_helper_integral(double x) {
    int position = get_index(x);
    if (position < v0plus_helper_values.size()) {
      return v0plus_helper_values[position];
    }
    throw std::out_of_range("v0plus helper cache miss");
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

struct increment {
  double next_point, delta_t;
  int change_grid_idx;
};
struct cellJob {
  std::unique_ptr<std::condition_variable> cv;
  std::unique_ptr<std::mutex> m;
  bool done;
};
template <typename DT> class double_vec {
public:
  void insert(int idx, DT x) {
    if ((idx >= 0) && (idx == right.size())) {
      right.push_back(std::move(x));
    } else if ((idx < 0) && (-1 - idx == left.size())) {
      left.push_back(std::move(x));
    } else if (!(((idx >= 0) && (idx < right.size())) ||
                 ((idx < 0) && (-1 - idx < left.size())))) {
      throw std::out_of_range("past end");
    }
  }
  DT &at(int idx) {
    if (idx >= 0) {
      return right.at(idx);
    } else {
      return left.at(-1 - idx);
    }
  }
  bool contains(int idx) {
    return (((idx >= 0) && (idx < right.size())) ||
            ((idx < 0) && (-1 - idx < left.size())));
  }

private:
  std::vector<DT> left;
  std::vector<DT> right;
};

class Simulator {
public:
  double start;
  Simulator(double x) { start = x; }
  std::vector<double> results;
  void simulate(double t, int rounds = ROUNDS,
                int print_interval = PRINT_INTERVAL) {
    std::vector<std::thread> threads;
    for (int idx = 0; idx < rounds; idx++) {
      threads.push_back(std::thread(&Simulator::run_sim, this, t, idx));
#ifdef _DEBUG
      std::cout << "spawned " << idx << std::endl;
#endif
    }
    for (int idx = 0; idx < rounds; idx++) {
      threads[idx].join();
#ifdef _DEBUG
      std::cout << "joined " << idx << std::endl;
#endif
    }
    std::ofstream output_file("res.csv");
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    std::copy(std::begin(results), std::end(results), output_iterator);
  }

private:
  double_vec<cellData> cell_cache;
  std::mutex cell_cache_mutex;
  std::mutex results_mutex;
  std::mutex cout_mutex;
  cellData get_data(double point, int grid_idx, int idx = 0) {
    std::unique_lock<std::mutex> lk(cout_mutex, std::defer_lock);

    if (cell_cache.contains(grid_idx)) { // cell in cache
#ifdef _DEBUG_VERBOSE
      lk.lock();
      std::cout << idx << " is at " << point << " cell cache hit" << std::endl;
      lk.unlock();
#endif
      return cell_cache.at(grid_idx);
    }
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " cell cache miss, computing"
              << std::endl;
    lk.unlock();
#endif
    cell lr = get_adjacent(point);
    CellDataCalculator calc(lr.left, lr.right, point);
    cellData out = calc.compute_cell_data();
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " finished job, writing"
              << std::endl;
    lk.unlock();
#endif
    std::unique_lock<std::mutex> cache_lock(cell_cache_mutex);
    cell_cache.insert(grid_idx, out);
    cache_lock.unlock();
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << ", writing and returning"
              << std::endl;
    lk.unlock();
#endif
    return out;
  }
  increment next_point(double point, std::mt19937 &rng, int grid_idx,
                       int idx = 0) {
    cell lr = get_adjacent(point);
    cellData point_data = get_data(point, grid_idx, idx);
    std::bernoulli_distribution d(point_data.prob_right);
    if (d(rng)) { // exit right
      return increment{lr.right, point_data.time_right, 1};
    } else { // exit left
      return increment{lr.left, point_data.time_left, -1};
    }
  }
  void run_sim(double t, int idx) {
    double cur = start;
    double t_cur = 0;
    std::unique_lock<std::mutex> lk(cout_mutex, std::defer_lock);
    std::random_device rd;
    std::mt19937 rng{rd()};

    int grid_idx = 0;

    while (t_cur < t) {
#ifdef _DEBUG_VERBOSE
      lk.lock();
      std::cout << idx << " is at " << cur << " at time " << t_cur << std::endl;
      lk.unlock();
#endif
      increment inc = next_point(cur, rng, grid_idx, idx);
      cur = inc.next_point;
      grid_idx += inc.change_grid_idx;
      t_cur += inc.delta_t;
    }

#ifdef _DEBUG
    lk.lock();
    std::cout << idx << " finished at " << cur << std::endl;
    lk.unlock();
#endif
    std::unique_lock<std::mutex> res_lock(results_mutex);
    results.push_back(cur);
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " pushed results " << std::endl;
    lk.unlock();
#endif
  }
};

int main() {
  Simulator sim(START);
  sim.simulate(TIME);
  return 0;
}
