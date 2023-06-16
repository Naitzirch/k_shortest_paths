# WIP
## What is this repository?
The code in this repositroy is part of my (Christian Martens's) bachelor final project.<br>
It aims to provide an implementation of Eppstein's k shortest paths algorithm (1997) using the networkx library.<br>
It was written for Patrick Emonts from the aQa group at Leiden University.<br>

## What can I expect from this code?
The code from the repository aims to implement the base algorithm of Eppstein's paper which should run in O(m + nlogn + klogk) (m edges, n vertices, k paths).<br>
The current implementation does not run within these time bounds.<br>
This implementation only accepts directed acyclic graphs (DAGs) with positive edge weights.<br>

## Author notes
Code in the epsnet folder was written by Patrick Emonts.
Code in test_graph and test_ising was also adapted from his code.
