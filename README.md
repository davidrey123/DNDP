# DNDP
Codes, instances and numerical results for the bilevel discrete network design problem (DNDP).

This repository contains implementation codes and benchmarking instances for the DNDP. The DNDP is a well-studied bilevel optimization problem in transportation. The goal of the DNDP is to identify the optimal set of candidate links (or projects) to be added to the network while accounting for users' reaction as governed by a traffic equilibrium model.

This repository complements the paper titled _Computational benchmarking of exact methods for the bilevel discrete network design problem_ by Rey, D. 2019.

__Implementation codes__ of exact methods for the linearized DNDP as discussed in the paper are provided in Python. These scripts require CPLEX's DOcplex Python API (https://pypi.org/project/docplex/) to solve MILP and LP problems, and the Pyomo module (http://www.pyomo.org/) along with IPOPT (https://www.coin-or.org/Ipopt/) to solve convex NLP problems. 

__Benchmark instances__ for the DNDP are provided for conceiving, testing and comparing solution methods. Each instance refers to a transport network publicly available at https://github.com/bstabler/TransportationNetworks in the TNTP format. Each instance contains the list of links of the original network (A_1) and the candidate links under consideration (A_2) along with their addition cost. The instances extend the TNTP format by including an extra column (right-most column) providing the cost of adding the corresponding link to the network (g_ij in the notation of the paper). 

__Detailed numerical results__ for all the numerical experiments conducted in the paper are provided in the supplementary material document. 
