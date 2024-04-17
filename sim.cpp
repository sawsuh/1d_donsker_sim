#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <fstream>
#include <future>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <stdexcept>
#include <unistd.h>
#include <unordered_map>

enum PlusMinus { plus, minus };

// INPUT
const int ROUNDS = 700;
const int PRINT_INTERVAL = 100;
const double INTEGRATION_INC = 0.0001;
const double START = 0;
const double TIME = 1;
double a(double x) { return 1; }
double rho(double x) { return 1; }
double b(double x) { return 0; }
struct cell {
  double left, right;
};
cell get_adjacent(double point) { return cell{point - 0.2, point + 0.2}; }

struct cellData {
  double time_left;
  double time_right;
  double prob_left;
  double prob_right;
};

void push_safe(std::vector<double> &vec, int idx, double x) {
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
#ifndef _SAFE
    throw std::out_of_range("psi cache miss");
#else
    double integral = 0;
    double y = left;
    while (y < x) {
      double y_next = y + integration_inc;
      int position_y = get_index(y_next);
      if (position_y < psi_values.size()) {
        integral = psi_values[position_y];
        y = y_next;
        continue;
      }
      integral +=
          integration_inc *
          (b(y) / (a(y) * rho(y)) + b(y_next) / (a(y_next) * rho(y_next))) / 2;
      push_safe(psi_values, position_y, integral);
      y = y_next;
    }
    push_safe(psi_values, position, integral);
    return integral * 2;
#endif
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
#ifndef _SAFE
    throw std::out_of_range("v0plus helper cache miss");
#else
    double integral = 0;
    double y = left;
    while (y < x) {
      double y_next = y + integration_inc;
      int position_y = get_index(y_next);
      if (position_y < v0plus_helper_values.size()) {
        integral = v0plus_helper_values[position_y];
        y = y_next;
        continue;
      }
      integral += integration_inc *
                  (exp(-psi(y)) / a(y) + exp(-psi(y_next)) / a(y_next)) / 2;
      push_safe(v0plus_helper_values, position_y, integral);
      y = y_next;
    }
    push_safe(v0plus_helper_values, position, integral);
    return integral;
#endif
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
      // threads[idx].join();
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
  std::unordered_map<int, cellData> cell_cache;
  std::mutex cell_cache_mutex;
  std::unordered_map<int, cellJob> current_cell_jobs;
  std::mutex current_cell_jobs_mutex;
  std::mutex results_mutex;
  std::mutex cout_mutex;
  cellData get_data(double point, int grid_idx, int idx = 0) {
    std::unique_lock<std::mutex> lk(cout_mutex, std::defer_lock);

    std::unique_lock<std::mutex> cache_lock(cell_cache_mutex);
    if (cell_cache.find(grid_idx) != cell_cache.end()) { // cell in cache
#ifdef _DEBUG_VERBOSE
      lk.lock();
      std::cout << idx << " is at " << point << " cell cache hit" << std::endl;
      lk.unlock();
#endif
      cache_lock.unlock();
      return cell_cache.at(grid_idx);
    }
    cache_lock.unlock();
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " cell cache miss" << std::endl;
    lk.unlock();
#endif
    std::unique_lock<std::mutex> jobs_lock(current_cell_jobs_mutex);
    if (current_cell_jobs.find(grid_idx) != current_cell_jobs.end()) {
      jobs_lock.unlock();
      // cell not in cache but job running
      std::unique_lock<std::mutex> job_wait_lock(
          *current_cell_jobs.at(grid_idx).m);
#ifdef _DEBUG_VERBOSE
      lk.lock();
      std::cout << idx << " is at " << point << " job ongoing, waiting"
                << std::endl;
      lk.unlock();
#endif
      while (!current_cell_jobs.at(grid_idx).done) {
        current_cell_jobs.at(grid_idx).cv->wait(job_wait_lock);
      }
#ifdef _DEBUG_VERBOSE
      lk.lock();
      std::cout << idx << " is at " << point << " waited job done" << std::endl;
      lk.unlock();
#endif
      cache_lock.lock();
      return cell_cache.at(grid_idx);
      cache_lock.unlock();
    }
    // cell not in cache and no job
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " no job ongoing, doing work"
              << std::endl;
    lk.unlock();
#endif

    // place job in joblist
    current_cell_jobs.insert(
        {grid_idx,
         cellJob{std::unique_ptr<std::condition_variable>(
                     new std::condition_variable),
                 std::unique_ptr<std::mutex>(new std::mutex), false}});
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " placed job in joblist"
              << std::endl;
    lk.unlock();
#endif
    jobs_lock.unlock();
    // job in joblist

    cell lr = get_adjacent(point);
    CellDataCalculator calc(lr.left, lr.right, point);
    cellData out = calc.compute_cell_data();
#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " finished job, writing"
              << std::endl;
    lk.unlock();
#endif
    cache_lock.lock();
    cell_cache.insert({grid_idx, out});
    cache_lock.unlock();

#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " written, notifying"
              << std::endl;
    lk.unlock();
#endif
    jobs_lock.lock();
    std::unique_lock<std::mutex> cur_job_lock(
        *current_cell_jobs.at(grid_idx).m);
    current_cell_jobs.at(grid_idx).done = true;
    current_cell_jobs.at(grid_idx).cv->notify_all();
    cur_job_lock.unlock();
    jobs_lock.unlock();

#ifdef _DEBUG_VERBOSE
    lk.lock();
    std::cout << idx << " is at " << point << " notified and erased"
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
