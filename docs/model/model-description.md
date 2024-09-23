
# Model Description

## Tree Structure

In the following, a summary of the model is provided. More detail is found in \cite Lignell_2024. The structure is based on a binary tree consisting of a hierarchy of nodes (represented by circles). At the base (bottom) of the tree, there are fluid parcels (represented by squares). The tree has five levels (indexed from 0 to 4). The top node, or apex, divides into two subtrees, each starting from a node one level below. Each of these nodes divides further into smaller subtrees. This pattern continues down to the base, where fluid properties are stored in parcels, and additional length and time scale information is associated with the nodes at each level.

For an N-level binary tree, the number of fluid parcels is \f$2^{N-1}\f$, and the length scale at each level is reduced by a factor of \f$A\f$ from the previous level, given by \f$L_{i+1} = L_iA\f$, which can be expressed as:

$$
L_i = L_0A^{i}
$$

where \f$A < 1\f$. This shows the multiplicative reduction in length scales as the tree progresses downward. In this model, the volume of each subtree is half that of its parent level. For three-dimensional space, we set \f$A = 1/2^{1/3}\f$, which maintains the volume halving along with the corresponding length scale reduction. In the general case, with \f$D\f$ dimensions, \f$A\f$ is calculated as \f$1/2^{1/D}\f$. In our study, we focus on a one-dimensional system, where \f$D = 1\f$, making \f$A = 1/2\f$. This leads to a significant reduction in length scale for each level, compared to higher \f$A\f$ values.

<div style="text-align: center;">
    <img src="hips_tree.png" style="width: 500px;">
</div>
From a mixing perspective, turbulent advection involves a gradual reduction in the size of flow structures through a cascading process that operates across different scales. According to Kolmogorov's second similarity hypothesis, within the inertial range of the flow, the mean kinetic energy dissipation rate, \f$\epsilon\f$, remains constant regardless of the length scale. This leads to the relation \f$\epsilon \sim l^2/\tau^3\f$, which gives \f$\tau \sim l^{2/3}\f$, and can be expressed as:

$$
\tau_i = \tau_0\left(\frac{L_i}{L_0}\right)^{2/3} = \tau_0A^{2i/3}.
$$

This relationship implies a weak coupling between different length scales, indicating that the breakdown of flow structures occurs gradually through small steps, a behavior observed empirically.

In this context, the scaling of velocity, \f$v \sim l^{1/3}\f$, suggests that larger fluid parcels exhibit higher velocities compared to smaller parcels. Consequently, the motions at larger scales influence smaller scales.

The Reynolds number, \f$Re\f$, is scaled based on Kolmogorov's first similarity hypothesis. It is related to the largest and smallest length scales in the inertial range, given by:

$$
Re = A^{-\frac{4}{3}i}.
$$

In the HiPS model, \f$Re\f$ is specified according to the number of levels in the system.

Turbulent mixing in HiPS is modeled by rearranging pairs of fluid parcels. The binary tree structure facilitates and defines this process. An "eddy event" describes the rearrangement, which is defined as follows:

1. A node in the tree is selected from levels 0 to N-3 (the base node).
2. Two random parcels (or nodes) are selected two levels down: one from the right branch and one from the left branch (grandchild parcels).
3. The subtrees beneath these two grandchild parcels are swapped.

Two types of swaps can occur. In the first case, the base node is selected at level 2. The grey-checked parcels labeled "a" and "d" are selected randomly and then swapped, changing their pairing from \f$(a, b)\rightarrow(a,d)\f$ and \f$(c, d)\rightarrow (c,a)\f$. In the second case, the swap happens at level 1, where the orange-hashed parcels \f$i\f$ and \f$j\f$ are swapped with parcels \f$m\f$ and \f$n\f$. 

Only swaps like the first example influence mixing directly, as they alter the parcel pairings. However, the other swaps still play a role by reducing the proximity-based variation in parcel properties, eventually enabling molecular mixing. At level \f$N-3\f$, parcel mixing corresponds to micromixing, while mixing at higher levels can be considered macromixing, affecting length scales without changing fine-scale fluid properties.


## Eddy Selection

In the HiPS tree, the eddy rate \f$\lambda_i\f$ at level \f$i\f$ is determined as \f$1/\tau_i\f$ for each node. Collectively, for all nodes at level \f$i\f$, the eddy rate is calculated as:

$$
\lambda_i = \frac{2^i}{\tau_i}.
$$

The occurrence times of eddy events are sampled from an exponential distribution, which represents a Poisson process with a mean rate \f$\Lambda\f$. The probability distribution is given by:

$$
p(\Delta t) = \Lambda e^{-\Lambda\Delta t}.
$$

To generate eddy occurrence times, we can sample from this distribution using the formula:

$$
\Delta t = -\frac{\ln(P_r)}{\Lambda},
$$

where \f$P_r\f$ denotes a uniform random variate in the range \([0,1]\).

To determine the tree level for a sampled eddy event, we follow a process that captures the full range of scalar fluctuations. These range from the inertial scales, \f$l \ge l^*\f$, denoted as \f$I\f$, to the viscous scales, \f$l \le l^*\f$, denoted as \f$V\f$. Here, \f$l^*\f$ is the characteristic length scale that marks the transition between the inertial and viscous regimes, and \f$l\f$ is any arbitrary length scale. 

The total rate is \f$\Lambda = \Lambda_I + \Lambda_V\f$. A uniform random variate \f$P_r\f$ is used to select the region. If \f$P_r \le \Lambda_I/\Lambda\f$, region \f$I\f$ is selected; otherwise, region \f$V\f$ is chosen, and then a specific level within the selected region is chosen.

Within the inertial range, the probability of an eddy event occurring at level \f$i\f$ is given by:

$$
p(i) = \frac{\lambda_i}{\Lambda_I},
$$

where \f$\Lambda_I = \sum_{i=0}^{i^*} \lambda_i\f$. The cumulative distribution function (CDF) is defined as:

$$
P(i) = \frac{1-(2A^{-2/3})^{i+1}}{1-(2A^{-2/3})^{i^*+1}}.
$$

The level \f$i\f$ can be sampled using:

$$
i = \left\lceil \frac{\log_2( 1-P_r(1-(2A^{-2/3})^{i^*+1}) )}{1-\frac{2}{3}\log_2A} - 1 \right\rceil.
$$

Within the viscous range, the probability of an eddy event occurring at level \f$i\f$ is given by:

$$
p(i) = \frac{\lambda_i}{\Lambda_V},
$$

where \f$\Lambda_V = \sum_{i=i^* +1}^{N-3}\lambda_i\f$. The CDF is expressed as:

$$
P(i) = \frac{2^{i+1}- 2^{i^* +1}}{2^{N-2}- 2^{i^* +1}}.
$$

The level \f$i\f$ can be sampled using a uniform random variate \f$P_r\f$ with the following formula:

$$
i = \left\lceil\log_2\left(P_r(2^{N-2}-2^{i^* +1})+2^{i^* +1}\right) -1 \right\rceil.
$$


##  Schmidt Dependence
The Schmidt number (Sc) represents the ratio of momentum diffusivity to species diffusivity. In our context, we extend this definition to apply to arbitrary scalars, not limited to chemical species. The variable Sc formulation in the HiPS model varies depending on whether Sc is greater than 1 or less than 1. This discrepancy arises from differences in how mixing is handled in the HiPS model compared to diffusion-based simulation approaches that typically model transport processes in physical coordinates. We first present the variable Sc model regarding scalars aligned with tree levels in length and timescales before extending it to arbitrary scales.

The \f$l^*\f$ scale holds particular significance for a scalar with Sc equal to unity and can be directly compared to the Kolmogorov scale \f$\eta\f$, with the two scales being proportional. Additionally, \f$l^*_s\f$ represents the smallest length scale for a scalar with an arbitrary \f$Sc\f$, akin to the Batchelor scale \f$\eta_b\f$ (for \f$Sc\ge 1\f$) or the Obukhov-Corrsin scale \f$\eta_{oc}\f$ (for \f$Sc\le 1\f$), and is proportional to \f$\eta_b\f$ or \f$\eta_{oc}\f$, respectively. The HiPS Schmidt number is defined as:

$$
Sc=(l^*/l_s^*)^{p_s}.
$$

For \f$Sc\ge 1\f$, the choice of \f$p_s=2\f$ is grounded in the analogy to a physical flow, where \f$\tau_\eta=\tau_{\eta_b}\f$ in the viscous regime. Here, \f$Sc=\nu/D\f$, \f$\nu=\eta^2/\tau_\eta\f$, and \f$D=\eta_b^2/\tau_{\eta_b}\f$.

For \f$Sc\le 1\f$, the selection of \f$p_s=3/4\f$ is similarly based on an analogy to a physical flow, where \f$\tau_{\eta_b}=\tau_\eta(\eta_b/\eta)^{2/3}\f$ in the inertial range (using \f$\tau_i = \tau_0A^{2i/3}\f$).


### Discretized Schmidt Values

To establish the relationship between \f$Sc\f$ and the associated scales in the HiPS tree, we initially focus on \f$l^*_s\f$ values restricted to HiPS levels, and then generalize it to arbitrary \f$l^*_s\f$. In a general HiPS simulation, multiple scalars with different \f$Sc\f$ values can be considered.

The \f$Sc\f$ value is related to the tree levels \f$i^*\f$ and \f$i_s^*\f$ as follows:

$$
Sc = A^{p_s(i^*-i^*_s)} = A^{-p_s\Delta i}.
$$

For \f$A=1/2\f$ and \f$p_s=4/3\f$ for \f$Sc\le 1\f$, we get \f$Sc=4^{2\Delta i/3}\f$, resulting in \f$Sc \approx\f$ 1, 0.4, 0.16, 0.062, 0.025, for \f$\Delta i=\f$ 0, -1, -2, -3, -4, respectively. For \f$Sc\ge 1\f$, with \f$p_s=2\f$, \f$Sc=4^{\Delta i}\f$, leading to \f$Sc=\f$ 1, 4, 16, 64, 256, for \f$\Delta i\f$ = 0, 1, 2, 3, 4, respectively.


### Arbitary Schmid


The formulation presented above identifies \f$Sc\f$ values associated with integer levels of the HiPS tree. Arbitrary \f$Sc\f$ corresponds to scalars with \f$l^*_s\f$ between two HiPS tree levels, and therefore, \f$i^*_s\f$ may not be an integer. For a scalar with a given \f$Sc\f$, \f$i^*_s\f$ is computed as 

$$
i^*_s = i^* - \frac{\log Sc}{p_s \log A},
$$ 

where:

- For \f$Sc < 1\f$:

$$
i^*_s = i^* - \frac{3 \log Sc}{4 \log A},
$$ 

- For \f$Sc > 1\f$:

$$
i^*_s = i^* - \frac{\log Sc}{2 \log A}.
$$

Levels \f$i_-\f$ and \f$i_+\f$ are considered for the lower and upper levels bounding \f$i^*_s\f$. Eddy events that occur on levels at or above \f$i_+\f$ mix the scalar across the left and right subtrees emanating from the \f$i_+\f$-level eddy node. For a level \f$i_-\f$ eddy event, the scalar is mixed across the two subtrees of the \f$i_-\f$-level node with probability \f$p_-\f$, where

$$
p_- = \frac{i_+ - i^*_s}{i_+ - i_-} = i_+ - i^*_s.
$$

The second equality holds since \f$i_+ - i_- \f$ is always unity. This probability is linear in index space and takes a value of 1 when \f$i^*_s = i_-\f$ and 0 when \f$i^*_s = i_+\f$. Using Equations for \f$L_i\f$, \f$\tau_i\f$, and \f$\lambda_i\f$, \f$p_-\f$ can be written as 

$$
p_- = \frac{\log(l^*_s/l_+)}{\log(l_- / l_+)} = \frac{\log(\lambda^*_s/\lambda_+)}{\log(\lambda_-/\lambda_+)}.
$$

This form illustrates that while there is a linear interpolation of \f$p_-\f$ in index space, the corresponding interpolation between eddy lengths, times, or rates is logarithmic, consistent with the geometric progression of the scales with tree level.

