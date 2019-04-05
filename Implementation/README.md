__Implementation codes for solving the bilevel discrete network design problem (DNDP)__

The file DNDP-Benchmark.py contains Python code for implementing three solution methods for the DNDP. These three methods are based on the following papers:

__SOBB__: Farvaresh, H. and Sepehri, M.M., 2013. A branch and bound algorithm for bi-level discrete network design problem. Networks and Spatial Economics, 13(1), pp.67-106.

__SOIC__: Wang, S., Meng, Q. and Yang, H., 2013. Global optimization methods for the discrete network design problem. Transportation Research Part B: Methodological, 50, pp.42-60.

__MKKT__: Fontaine, P. and Minner, S., 2014. Benders decomposition for discrete–continuous linear bilevel problems with application to traffic network design. Transportation Research Part B: Methodological, 70, pp.163-172.

All three solution methods, SOBB, SOIC and MKKT, are implemented using the piece-wise linear approximation of link delay functions proposed by Farvaresh, H. and Sepehri, M.M., 2011. A single-level mixed integer linear formulation for a bi-level discrete network design problem. Transportation Research Part E: Logistics and Transportation Review, 47(5), pp.623-640.
In all cases, the optimization problems solved are MILPs and LPs (traffic assignment problems within SOBB and SOIC are solved in their linearized form).

 All MILPs and LPs are solved using CPLEX 12.8 MIP solver. Convex traffic assignment problems (for verification) are solved using the Pyomo module and IPOPT solver.
