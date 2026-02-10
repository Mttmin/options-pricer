# Overall To-do list

- Convertible pricing using more than BS
- ✅ American puts (PenaltySolver working, 1-2% accuracy vs binomial)
- American calls
- Asian Call/puts
- All around multi-threaded Monte Carlo -> adapt to American and Asian options
- Auto-callables pricing
- Greeks calculation for exotics based on Monte-carlo
- API for live info of option prices and testing compared to market data on certain SPY options
- Web scrapping/api access for: Automatic Stock data info,Vol calculation, and Risk-free rate pull

## American Options PenaltySolver - Performance Optimization

- [x] Phase 1: Fixed critical bugs and added diagnostics (COMPLETE)
  - Fixed drift term sign error (was -r*S, now +rS)
  - Extended spatial grid (10σ width, 0.1% minimum)
  - Added solver diagnostics
  - Achieves 1.33% accuracy at 800x400 grid

- [x] Phase 2: Optimizations implemented (COMPLETE)
  - ✅ Precomputed spatial coefficients (parallel with rayon)
  - ✅ Vectorized penalty diagonal construction
  - ✅ Parallel convergence checking
  - ✅ Cached base matrices across iterations
  - Note: Full sparse matrix solver deferred (requires external library like MUMPS)

  Performance results (release mode):
  - 100×50:  21ms (1.68% error)
  - 200×100: 127ms (4.90% error)
  - 400×200: 1.4s (3.05% error)
  - 800×400: 20.8s (1.33% error)

- [ ] Phase 3: Better test cases and validation
  - Match paper's Table 8.1 results
  - Broader parameter tests
  - Test American calls

- [ ] Future optimizations (when needed for N > 1000):
  - Iterative solvers (BiCGSTAB/GMRES) via external library
  - GPU acceleration for large grids
  - Multigrid methods
