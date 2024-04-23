#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <omp.h>
#include <random>
#include <set>
#include <stdexcept>
#include <thread>
#include <unistd.h>
#include <unordered_map>

enum PlusMinus { plus, minus };

// INPUT

// Number of simulations
const int ROUNDS = 10000;
// Interval to use for numerical integration
const double INTEGRATION_INC = 0.00001;
// Starting point in grid
const double START = 0;
// Time at which we want to sample X_t
const double TIME = 1;
// a
double a(double x) { return 1; }
// rho
double rho(double x) {
  if (x < 0) {
    return 1;
  }
  return 5;
}
// b
double b(double x) { return 0; }
// =============================

// cell: contains values of left and right neighbours
struct cell {
  double left, right;
};
// GRID SPECIFICATION
// args: point in grid
// returns: two adjacent points
cell get_adjacent(double point) { return cell{point - 0.01, point + 0.01}; }
// =============================

// data for a cell:
// contains left and right exit probabilities and
//  left and right exit times conditional on left and right exit probabilities
struct cellData {
  double time_left;
  double time_right;
  double prob_left;
  double prob_right;
};

// vector that we use to cache data for our integral
template <typename DT> void push_safe(std::vector<DT> &vec, int idx, DT x) {
  if (vec.size() != idx) {
    throw std::out_of_range("desired location not end of vector");
  }
  vec.push_back(x);
}

// handles computation of integrals for a single cell
class CellDataCalculator {
public:
  // left boundary
  double left;
  // right boundary
  double right;
  // the point in the cell which our process starts at
  double x;
  CellDataCalculator(double l, double r, double y) {
    left = l;
    right = r;
    x = y;
    integration_inc = INTEGRATION_INC;
  }
  // calculate data (exit times and probabilities)
  cellData compute_cell_data() {
    // precompute table of psi values
    gen_psi_table();
    // precompute table of s values
    gen_s_table();
    // we can now compute exit probabilities
    double prob_left = v0minus(x);
    double prob_right = v0plus(x);
    // we can now compute exit times
    double time_left_ind = v1minus(x);
    double time_right_ind = v1plus(x);
    return cellData{time_left_ind / prob_left, time_right_ind / prob_right,
                    prob_left, prob_right};
  }

private:
  // table of values of psi
  std::vector<double> psi_values;
  // table of values of s
  std::vector<double> v0plus_helper_values;
  // increment used for numerical integration
  double integration_inc;
  // arg: point
  // returns: index in table corresponding to the point (starts at 0)
  // (for table access)
  int get_index(double x) { return round((x - left) / integration_inc); }
  // progressively precompute values of psi in one pass
  void gen_psi_table() {
    // reserve vector (we know the necessary length)
    psi_values.reserve(get_index(right));
    // integral from left to left is 0
    psi_values.push_back(0);
    double integral = 0;
    double y = left;
    double y_next;
    // integrate step by step and store each step as the integral up to that
    // point
    while (y < right) {
      y_next = y + integration_inc;
      // use trapezoid rule to integrate
      integral += integration_inc * (b(y) / (a(y) * rho(y)) +
                                     b(y_next) / (a(y_next) * rho(y_next)));
      // store integral
      push_safe(psi_values, get_index(y_next), integral);
      y = y_next;
    }
  }
  // all the necessary values should be stored in cache, throw an error if not
  double psi(double x) {
    int position = get_index(x);
    if (position < psi_values.size()) {
      return psi_values[position];
    }
    throw std::out_of_range("psi cache miss");
  }
  // precompute values of s, same strategy as psi
  // (do integral in one pass, storing each step as the integral to that point)
  void gen_s_table() {
    // we know the full size
    v0plus_helper_values.reserve(get_index(right));
    // integral left to left is 0
    v0plus_helper_values.push_back(0);
    double integral = 0;
    double y = left;
    double y_next;
    // for each step
    while (y < right) {
      y_next = y + integration_inc;
      // compute integral to next point by trapezoid rule
      integral += integration_inc *
                  (exp(-psi(y)) / a(y) + exp(-psi(y_next)) / a(y_next)) / 2;
      // store
      push_safe(v0plus_helper_values, get_index(y_next), integral);
      y = y_next;
    }
  }
  // all values should be precomputed, if not throw an error
  double s(double x) {
    int position = get_index(x);
    if (position < v0plus_helper_values.size()) {
      return v0plus_helper_values[position];
    }
    throw std::out_of_range("v0plus helper cache miss");
  }
  // Now that we have s, we can compute v0plus
  double v0plus(double x) { return s(x) / s(right); }
  // v0minus = 1-v0plus
  double v0minus(double x) { return 1 - v0plus(x); }
  // We can also express G using s
  double G(double x, double y) {
    if (x <= y) {
      return 2 * (s(x) - s(left)) * (s(right) - s(y)) / (s(right) - s(left));
    } else {
      return 2 * (s(y) - s(left)) * (s(right) - s(x)) / (s(right) - s(left));
    }
  }
  // We can use G together with psi and v0plus to compute v1plus
  // We only do this once so we only need one pass
  // We use the trapezoid rule
  double v1plus(double x) {
    double integral = 0;
    for (int y_idx = 0; y_idx < get_index(right); y_idx++) {
      double y = y_idx * integration_inc + left;
      double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (G(x, y) * v0plus(y) * exp(psi(y)) / rho(y) +
           G(x, y_next) * v0plus(y_next) * exp(psi(y_next)) / rho(y_next)) /
          2;
    }
    return integral;
  }
  // We can use G together with psi and v0minus to compute v1minus
  // We only do this once so we only need one pass
  // We use the trapezoid rule
  double v1minus(double x) {
    double integral = 0;
    for (int y_idx = 0; y_idx < get_index(right); y_idx++) {
      double y = y_idx * integration_inc + left;
      double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (G(x, y) * v0minus(y) * exp(psi(y)) / rho(y) +
           G(x, y_next) * v0minus(y_next) * exp(psi(y_next)) / rho(y_next)) /
          2;
    }
    return integral;
  }
};

// stores an increment (when simulating random walk)
// we need to know:
//  the next point we go to
//  how long it takes us to get there
//  we also track where our index in the grid changes
struct increment {
  double next_point, delta_t;
  int change_grid_idx;
};
// this is a structure for storing data we compute correponding to points in the
// grid
class cell_cache {
public:
  std::atomic<int> left;
  std::atomic<int> right;
  std::mutex left_m;
  std::mutex right_m;
  std::list<cellData> container;
  cell_cache(double start) {
    left = 0;
    right = 0;
    // get neighbours
    cell lr = get_adjacent(start);
    // Compute values
    CellDataCalculator calc(lr.left, lr.right, start);
    container.push_back(calc.compute_cell_data());
    it = container.begin();
  }
  std::list<cellData>::iterator get_iter() { return it; }
  // void insert(int idx, cellData x) {}
  // cellData at(int idx) {}
  // bool contains(int idx) {}
private:
  std::list<cellData>::iterator it;
};
class cell_cache_walker {
public:
  cell_cache_walker(std::shared_ptr<cell_cache> &c) : cache(c) {
    it = c->get_iter();
    idx = 0;
    started = true;
  }
  bool contains(int goal) {
    if (started && (abs(goal - idx) <= 1)) {
      started = false;
      return ((goal >= (*cache).left) && (goal <= (*cache).right));
    } else if ((!started) && (abs(goal - idx) == 1)) {
      return ((goal >= (*cache).left) && (goal <= (*cache).right));
    }
    throw std::out_of_range("checking too far away");
  }
  cellData at(int goal) {
    if (((goal < (*cache).left) || (goal > (*cache).right)) ||
        (abs(goal - idx) > 1)) {
      throw std::out_of_range("checking more than a step away or out of range");
    }
    int change = goal - idx;
    idx += change;
    if (change == 1) {
      return *(++it);
    } else if (change == -1) {
      return *(--it);
    } else if (change == 0) {
      return *it;
    }
    throw std::out_of_range("checking wrong place");
  }
  bool insert(int goal, cellData x) {
    bool wrote = false;
    if (goal == idx - 1) {
      std::unique_lock<std::mutex> lk((*cache).left_m, std::defer_lock);
      if (lk.try_lock() && (goal == (*cache).left - 1)) {
        (*cache).container.push_front(x);
        (*cache).left -= 1;
        wrote = true;
      }
      idx = idx - 1;
      --it;
    } else if (goal == idx + 1) {
      std::unique_lock<std::mutex> lk((*cache).right_m, std::defer_lock);
      if (lk.try_lock() && (goal == (*cache).right + 1)) {
        (*cache).container.push_back(x);
        (*cache).right += 1;
        wrote = true;
      }
      idx = idx + 1;
      ++it;
    } else {
      throw std::out_of_range("inserting more than 1 away");
    }
    return wrote;
  }

private:
  std::shared_ptr<cell_cache> &cache;
  std::list<cellData>::iterator it;
  int idx;
  bool started;
};

// Handles simulating the random walks
class Simulator {
public:
  // Start point
  double start;
  Simulator(double x) : cache_ptr(std::make_shared<cell_cache>(x)) {
    start = x;
  }
  // Holds results
  std::vector<double> results;
  // Runs simulation
  void simulate(double t, int rounds = ROUNDS) {
    std::vector<std::thread> threads;
    for (int idx = 0; idx < rounds; idx++) {
      threads.push_back(std::thread(&Simulator::run_sim, this, t, idx));
    }
    for (int idx = 0; idx < rounds; idx++) {
      threads[idx].join();
    }
    // Write results to "res.csv"
    std::ofstream output_file("res.csv");
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    std::copy(std::begin(results), std::end(results), output_iterator);
  }

private:
  std::mutex cout_m;
  // Holds computed data for cells we have visited
  // Uses double vector because we visit cells from start point outwards
  // And can just count the jumps
  std::shared_ptr<cell_cache> cache_ptr;
  cellData get_data(cell_cache_walker &walker, double point, int grid_idx,
                    int idx) {
    std::unique_lock<std::mutex> lk(cout_m, std::defer_lock);
// If the point is in our cache
// (we already visited it and have data)
#ifdef _DEBUG
    lk.lock();
    std::cout << idx << " checking for " << grid_idx << std::endl;
    lk.unlock();
#endif
    if (walker.contains(grid_idx)) {
#ifdef _DEBUG
      lk.lock();
      std::cout << idx << " found " << grid_idx << std::endl;
      lk.unlock();
#endif
      return walker.at(grid_idx);
    }
#ifdef _DEBUG
    lk.lock();
    std::cout << idx << " couldnt find " << grid_idx << std::endl;
    lk.unlock();
#endif
    // If the point is not in our cache
    // Compute the necessary data

    // get neighbours
    cell lr = get_adjacent(point);
    // Compute values
    CellDataCalculator calc(lr.left, lr.right, point);
    cellData out = calc.compute_cell_data();
// write computed data to cache
#ifdef _DEBUG

    lk.lock();
    std::cout << idx << " trying to write " << grid_idx << std::endl;
    lk.unlock();
#endif
    bool written = walker.insert(grid_idx, out);
#ifdef _DEBUG
    if (written) {
      lk.lock();
      std::cout << idx << " wrote " << grid_idx << std::endl;
      lk.unlock();
    } else {
      lk.lock();
      std::cout << idx << " couldnt write " << grid_idx << std::endl;
      lk.unlock();
    }
#endif

    // return computed value
    return out;
  }
  // take a step from a point
  increment next_point(cell_cache_walker &walker, double point,
                       std::mt19937 &rng, int grid_idx, int idx = 0) {
    // get neighbours
    cell lr = get_adjacent(point);
    // get data
    cellData point_data = get_data(walker, point, grid_idx, idx);
    // simulate step
    std::bernoulli_distribution d(point_data.prob_right);
    if (d(rng)) {
      // exit right
      return increment{lr.right, point_data.time_right, 1};
    } else {
      // exit left
      return increment{lr.left, point_data.time_left, -1};
    }
  }
  // run 1 simulation
  // takes start point
  // idx represents iteration number (for logging)
  void run_sim(double t, int idx) {
    double cur = start;
    double t_cur = 0;
    // random device for random number generation
    std::random_device rd;
    std::mt19937 rng{rd()};
    cell_cache_walker walker(cache_ptr);

    int grid_idx = 0;

    // run algorithm
    while (t_cur < t) {
      // get next point
      increment inc = next_point(walker, cur, rng, grid_idx, idx);
      // update position in grid and increment time
      cur = inc.next_point;
      grid_idx += inc.change_grid_idx;
      t_cur += inc.delta_t;
    }
#ifdef _DEBUG
    std::cout << idx << " finished at " << cur << "\n";
#endif
    // write result
    results.push_back(cur);
  }
};

int main() {
  Simulator sim(START);
  sim.simulate(TIME);
  return 0;
}
