# Testing repset's optima

Finding an optimal (minimizing the sum of percent identities) representative set of sequences is NP-hard. It can be solved exactly with the following MILP formulation:

$\min\limits_{x, y} \; \sum\limits_{i}\sum\limits_{j > i} I_{ij} y_{ij}$ <br>
s.t. <br>
<br>
$x_i + x_j - 2 y_{ij} \leq 1 \;\; (1)$ <br>
$x_i + x_j - 2 y_{ij} \geq 0 \;\; (2)$ <br>
$\sum\limits_{i} x_i = K \;\; (3)$ <br>
<br>
with <br>
$x \in \{0, 1\}^n$ <br>
$y \in \{0, 1\}^{\frac{n(n-1)}{2}}$ <br>

where $I_{k_{ij}}$ corresponds to the percent identity value between sequence $i$ and $j$.

While the MILP approach is guaranteed to reach a global optimum, running a MILP becomes impractical when there are more than a few hundred query sequences.

Instead of the previous MILP, repset relies on a submoldular optimization approach, which is scalable to thousands of query sequences. However, it's solution does not necessarily correspond to the global optimum of the previous MILP.

Here, we will test repset's optimal solution agains that of th MILP formulation to assess their average distance.