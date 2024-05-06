#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <random>

#pragma omp declare reduction(                                                 \
        vec_concat : std::vector<long double> : omp_out.insert(                \
                omp_out.end(), omp_in.begin(), omp_in.end()))                  \
    initializer(omp_priv = omp_orig)

// represents left or right exit
enum PlusMinus { plus, minus };
// cell: contains values of left and right neighbours
struct cell {
  const long double left, right;
};

// INPUT

// Number of simulations
const int ROUNDS = 10000;
// Interval to use for numerical integration
const long double INTEGRATION_INC = 0.00001;
// default: Brownian motion
#if !(defined(_FIG1) || defined(_FIG2) || defined(_FIG3) || defined(_FIG4) ||  \
      defined(_FIG5))
// Starting point in grid
const long double START = 0;
// Time at which we want to sample X_t
const long double TIME = 1;
// a
long double a(const long double x) { return 1; }
// rho
long double rho(const long double x) { return 1; }
// b
long double b(const long double x) { return 0; }
// GRID SPECIFICATION
// args: point in grid
// returns: two adjacent points
cell get_adjacent(long double point) {
  return cell{point - 0.01, point + 0.01};
}
#endif

// defining params for each of the figs
#if (defined(_FIG1) || defined(_FIG2) || defined(_FIG3))
const long double START = 1;
long double a(long double x) {
  if (x < -1.5) {
    return 1;
  } else if (x < 0.5) {
    return 2;
  } else {
    return 1;
  }
}
long double rho(long double x) {
  if (x < -0.5) {
    return 1;
  } else if (x < 0.5) {
    return 0.5;
  } else {
    return 1;
  }
}
long double b(long double x) { return 0; }
cell get_adjacent(long double point) {
  return cell{point - 0.02, point + 0.02};
}
#endif
#if defined(_FIG1)
const long double TIME = 0.5;
#elif defined(_FIG2)
const long double TIME = 1;
#elif defined(_FIG3)
const long double TIME = 1.5;
#endif
#if (defined(_FIG4) || defined(_FIG5))
const long double START = 0;
const long double TIME = 1;
long double b(long double x) { return 0; }
cell get_adjacent(long double point) {
  return cell{point - 0.05, point + 0.05};
}
#endif

#if defined(_FIG4)
long double a(long double x) { return (x < 0) ? 1 : 5; }
long double rho(long double x) { return 1; }
#elif defined(_FIG5)
long double a(long double x) { return 1; }
long double rho(long double x) { return (x < 0) ? 1 : 5; }
#endif

// =============================

// data for a cell:
// contains left and right exit probabilities and
//  left and right exit times conditional on left and right exit probabilities
struct cellData {
  const long double time_left;
  const long double time_right;
  const long double prob_left;
  const long double prob_right;
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
  // boundary and start point
  const long double left, right, x;
  CellDataCalculator(long double l, long double r, long double y)
      : left(l), right(r), x(y), integration_inc(INTEGRATION_INC) {}
  // calculate data (exit times and probabilities)
  cellData compute_cell_data() {
    // precompute table of psi values
    gen_psi_table();
    // precompute table of s values
    gen_s_table();
    // we can now compute exit probabilities
    const long double prob_left = v0minus(x);
    const long double prob_right = v0plus(x);
    // we can now compute exit times
    const long double time_left_ind = v1minus(x);
    const long double time_right_ind = v1plus(x);
    return cellData{time_left_ind / prob_left, time_right_ind / prob_right,
                    prob_left, prob_right};
  }

private:
  // tables of values for psi and s
  std::vector<long double> psi_values, s_values;
  // increment used for numerical integration
  const long double integration_inc;
  // arg: point
  // returns: index in table corresponding to the point (starts at 0)
  // (for table access)
  int get_index(long double x) const {
    return round((x - left) / integration_inc);
  }
  // progressively precompute values of psi in one pass
  void gen_psi_table() {
    // reserve vector (we know the necessary length)
    psi_values.reserve(get_index(right) + 1);
    // integral from left to left is 0
    psi_values.push_back(0);
    long double integral = 0;
    long double y = left;
    long double y_next = left;
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
  long double psi(long double x) const { return psi_values.at(get_index(x)); }
  // precompute values of s, same strategy as psi
  // (do integral in one pass, storing each step as the integral to that point)
  void gen_s_table() {
    // we know the full size
    s_values.reserve(get_index(right) + 1);
    // integral left to left is 0
    s_values.push_back(0);
    long double integral = 0;
    long double y = left;
    long double y_next = left;
    // for each step
    while (y < right) {
      y_next = y + integration_inc;
      // compute integral to next point by trapezoid rule
      integral += integration_inc *
                  (exp(-psi(y)) / a(y) + exp(-psi(y_next)) / a(y_next)) / 2;
      // store
      push_safe(s_values, get_index(y_next), integral);
      y = y_next;
    }
  }
  // all values should be precomputed, if not throw an error
  long double s(long double x) const { return s_values.at(get_index(x)); }
  // Now that we have s, we can compute v0plus
  long double v0plus(long double x) const { return s(x) / s(right); }
  // v0minus = 1-v0plus
  long double v0minus(long double x) const { return 1 - v0plus(x); }
  // We can also express G using s
  long double G(long double x, long double y) const {
    if (x <= y) {
      return 2 * (s(x) - s(left)) * (s(right) - s(y)) / (s(right) - s(left));
    } else {
      return 2 * (s(y) - s(left)) * (s(right) - s(x)) / (s(right) - s(left));
    }
  }
  // We can use G together with psi and v0plus to compute v1plus
  // We only do this once so we don't need to cache
  // We use the trapezoid rule
  long double v1plus(long double x) const {
    long double integral = 0;
    for (int y_idx = 0; y_idx < get_index(right); y_idx++) {
      long double y = y_idx * integration_inc + left;
      long double y_next = y + integration_inc;
      integral +=
          integration_inc *
          (G(x, y) * v0plus(y) * exp(psi(y)) / rho(y) +
           G(x, y_next) * v0plus(y_next) * exp(psi(y_next)) / rho(y_next)) /
          2;
    }
    return integral;
  }
  // We can use G together with psi and v0minus to compute v1minus
  // We only do this once so we don't need to cache
  // We use the trapezoid rule
  long double v1minus(long double x) const {
    long double integral = 0;
    for (int y_idx = 0; y_idx < get_index(right); y_idx++) {
      long double y = y_idx * integration_inc + left;
      long double y_next = y + integration_inc;
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
  const long double next_point, delta_t;
  const int change_grid_idx;
};
// this is a structure for storing data we compute correponding to points in the
// grid
class cell_cache {
public:
  // actually holds the data
  std::list<cellData> container;
  // initialise with data for the start point
  cell_cache(long double start) : left(0), right(0) {
    // get neighbours
    cell lr = get_adjacent(start);
    // Compute values
    CellDataCalculator calc(lr.left, lr.right, start);
    // store
    container.push_back(calc.compute_cell_data());
    // save iterator to start
    it = container.begin();
  }
  // return iterator to start
  std::list<cellData>::iterator get_iter() const { return it; }
  // get left endpoint
  int get_left() const {
    int l;
    // atomic for thread safety
#pragma omp atomic read
    l = left;
    return l;
  }
  // get right endpoint
  int get_right() const {
    int r;
    // atomic for thread safety
#pragma omp atomic read
    r = right;
    return r;
  }
  // advance left endpoint
  void inc_left() {
    // atomic for thread safety
#pragma omp atomic update
    left--;
  }
  // advance right endpoint
  void inc_right() {
    // atomic for thread safety
#pragma omp atomic update
    right++;
  }

private:
  // persistent iterator to start point
  std::list<cellData>::iterator it;
  // endpoints
  int left, right;
};
// object that each thread uses to walk the cache
// effectively a wrapper around an iterator
class cell_cache_walker {
public:
  // uses a shared pointer to the cache
  cell_cache_walker(std::shared_ptr<cell_cache> c) : cache(c) {
    // start at start point
    it = cache->get_iter();
    idx = 0;
  }
  // check if cache holds an item
  // left and right are the ends of the cache
  // these are already atomic so no race condition here
  bool contains(int goal) const {
    return (cache->get_left() <= goal) && (goal <= cache->get_right());
  }
  // return an item that is 0 or 1 steps away
  // (and walk there)
  // assumes we already checked existence
  cellData at(int goal) {
    const int change = goal - idx;
    // increment internal position
    idx += change;
    // update iterator
    if (change == 1) {
      return *(++it);
    } else if (change == -1) {
      return *(--it);
    } else if (change == 0) {
      return *it;
    }
    throw std::out_of_range("indexing too far away");
  }
  // safely insert an item
  void insert(int goal, cellData x) {
    const int change = goal - idx;
    // left frontier
    if (change == -1) {
      // insert if it isn't there
      // critical to prevent races
#pragma omp critical(cache)
      {
        if (goal == cache->get_left() - 1) {
          cache->container.push_front(x);
          cache->inc_left();
        }
      }
      // walk this instance
      idx = idx - 1;
      it--;
    } else if (change == 1) {
      // right frontier
      // insert if not there
      // critical to prevent races
#pragma omp critical(cache)
      {
        if (goal == cache->get_right() + 1) {
          cache->container.push_back(x);
          cache->inc_right();
        }
      }
      // walk this instance
      idx = idx + 1;
      it++;
    } else {
      throw std::out_of_range("inserting not exactly 1 away");
    }
  }

private:
  // stores actual cache
  const std::shared_ptr<cell_cache> cache;
  // iterator that we walk around
  std::list<cellData>::iterator it;
  // tracks where in the grid our iterator is
  int idx;
};

// Handles simulating the random walks
class Simulator {
public:
  // Start point
  const long double start;
  // initialise cache
  Simulator(long double x)
      : start(x), cache_ptr(std::make_shared<cell_cache>(x)) {}
  // Runs simulation
  void simulate(long double t, int rounds = ROUNDS) const {
    // multithread this
    std::vector<long double> results;
#pragma omp parallel for reduction(vec_concat : results)
    for (int idx = 0; idx < rounds; idx++) {
      results.push_back(run_sim(t, idx));
    }
    // Write results to "res.csv"
    std::ofstream output_file("res.csv");
    std::ostream_iterator<long double> output_iterator(output_file, "\n");
    std::copy(std::begin(results), std::end(results), output_iterator);
#ifdef _DEBUG
    std::cout << "wrote " << results.size() << " results\n";
#endif
  }

private:
  // Holds computed data for cells we have visited
  const std::shared_ptr<cell_cache> cache_ptr;
  // computes values
  cellData get_data(cell_cache_walker &walker, long double point,
                    int grid_idx) const {
    // If the point is in our cache
    // (we already visited it and have data)
    if (walker.contains(grid_idx)) {
      return walker.at(grid_idx);
    }
    // If the point is not in our cache
    // Compute the necessary data

    // get neighbours
    const cell lr = get_adjacent(point);
    // Compute values
    CellDataCalculator calc(lr.left, lr.right, point);
    const cellData out = calc.compute_cell_data();
    // write computed data to cache
    walker.insert(grid_idx, out);

    // return computed value
    return out;
  }
  // take a step from a point
  increment next_point(cell_cache_walker &walker, long double point,
                       std::mt19937 &rng, int grid_idx) const {
    // get neighbours
    const cell lr = get_adjacent(point);
    // get data
    const cellData point_data = get_data(walker, point, grid_idx);
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
  // run 1 thread
  // takes start point
  // idx represents iteration number (for logging)
  long double run_sim(long double t, int idx) const {
#ifdef _DEBUG
#pragma omp critical(cout)
    std::cout << idx << " started\n";
#endif
    // initialise cache walker for this walk
    cell_cache_walker walker(cache_ptr);

    // initialise data
    long double cur = start;
    long double t_cur = 0;
    int grid_idx = 0;

    // random device for random number generation
    std::random_device rd;
    std::mt19937 rng{rd()};

    // run algorithm
    while (t_cur < t) {
      // get next point
      increment inc = next_point(walker, cur, rng, grid_idx);
      // update position in grid and increment time
      cur = inc.next_point;
      grid_idx += inc.change_grid_idx;
      t_cur += inc.delta_t;
    }
#ifdef _DEBUG
#pragma omp critical(cout)
    std::cout << idx << " finished at " << cur << "\n";
#endif
    return cur;
  }
};

int main() {
  Simulator sim(START);
  sim.simulate(TIME);
  return 0;
}
