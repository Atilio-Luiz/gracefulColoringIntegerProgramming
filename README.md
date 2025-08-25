# Graceful Coloring Integer Programming
Integer Programming for Graceful Coloring of Graphs

**Definition:** Let $G = (V(G),E(G))$ be simple graph $G$ and let $S = \{1,2,\ldots,k\}$ be a set of natural numbers called **colors**. A **k-graceful coloring** of $G$ is a function $f \colon V(G) \to \{1,2,\ldots,k\}$ such that the following two conditions are satisfied:

- every two adjacent vertices $u, v \in V(G)$ are assigned differente colors, that is, $f(u) \neq f(v)$; and
- when every edge $uv \in E(G)$ is assigned the induced edge label $|f(u)-f(v)|$, we obtain that any two adjacent edges have distinct induced edge labels.

The minimum whole number $k$ for which a graph $G$ admits a $k$-graceful coloring is called the **graceful chromatic number** of $G$ and is denoted by $\chi_g(G)$. The **Graceful Coloring Problem** consists in determining the exact value of the graceful chromatic number of an arbitrary graph. This problem is an NP-hard problem.

In this project, I propose an Integer Programming formulation for the graceful coloring problem and implement it in the framework provided by the Julia language, using CPLEX as the solver.