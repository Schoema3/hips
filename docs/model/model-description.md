
# Model Description

## Tree Structure

In the following, a summary of the model is provided. For more details, including a detailed analysis of the physical results of turbulent mixing, users can refer to \cite Lignell2024. Here, a schematic of the tree is shown.
<div style="text-align: center;">
    <img src="hips_tree.png" style="width: 500px;">
</div>

The tree consists of a hierarchy of nodes represented by circles in the figure. At the base (bottom) of the tree, there are fluid parcels represented by squares. The figure displays five levels (indexed 0 to 4). The top node (apex) divides into two subtrees, each starting from a node one level below the apex. Each of these nodes further divides into two smaller subtrees. This pattern continues down to the tree's base. Fluid properties are stored in parcels, while additional length and time scale information are linked to nodes at each level. For an N-level binary tree, the number of fluid parcels is \f$2^{N-1}\f$, and the length scale at each level is a factor \f$A\f$ of the length scale at the previous level, given by \f$L_{i+1} = L_iA\f$, 
which gives 
$$
L_i = L_0A^{i},
$$
\label{e:L}

where \f$A < 1\f$. This corresponds to successive multiplicative length scale reductions as the tree base approaches. In a binary tree model where parcels represent fluid volume, each sub-tree is half the volume of its parent tree level. If we set \f$A\f$ to \f$1/2^{1/3}\f$, it matches the volume halving with a proportionate reduction in length scale, suitable for a three-dimensional space. Generally, for a space with \f$D\f$ dimensions, \f$A\f$ is calculated as \f$1/2^{1/D}\f$. Our study focuses on a one-dimensional case, where \f$D=1\f$, making \f$A=1/2\f$. This choice leads to a more significant reduction in length scale for each level and parcel, compared to higher \f$A\f$ values. For additional details, see \cite Kerstein_2013.


From a mixing perspective, turbulent advection involves a gradual reduction in the size of flow structures through a scale-invariant cascading process that operates locally across various scales. Kolmogorov's second similarity hypothesis \cite Pope_2000 says that within the inertial range of the flow, the mean kinetic energy dissipation rate, \f$\epsilon\f$, remains constant regardless of the length scale, leading to \f$\epsilon \sim l^2/{\tau^3}\f$ and consequently 
\f$\tau \sim l^{2/3}\f$, 
which gives

$$
\tau_i = \tau_0\left(\frac{L_i}{L_0}\right)^{2/3} = \tau_0A^{2i/3}.
$$

This relationship suggests a weak coupling between different length scales, indicating that the breakdown of flow structures occurs gradually through a series of small steps, a behavior supported by empirical observations.

In this concept, the scaling of velocity, \f$v \sim l^{1/3}\f$, implies that larger fluid parcels at larger scales exhibit higher velocities \f$v\f$ compared to smaller parcels at smaller scales. Consequently, smaller scales are influenced and swept by the motions of larger scales.

The Reynolds number, represented as \f$Re\f$, is determined by Kolmogorov's first similarity hypothesis \cite Pope_2000. This hypothesis scales the Reynolds number relative to the largest and smallest length scales within the inertial range. Following this definition, we have:

$$
Re = A^{-\frac{4}{3}i}.
$$
Thus, in the HiPS model, \f$Re\f$ is specified according to a predetermined number of levels.

In HiPS, turbulent mixing is modeled by rearranging parcel pairs. The binary tree facilitates and defines the rearrangement process. The rearrangement, which is called an "eddy event," is defined as follows:

1. A given node of the tree is selected within accessible levels 0 to N-3. This is the base node.
2. A random node (or fluid parcel) two levels down on the right branch is selected, and another random node (or fluid parcel) is selected two levels down on the left branch. These are grandchild nodes/parcels.
3. The subtrees at and below the two grandchild parcels are swapped.

Two kinds of swapping processes are shown in the figure. The first is shown where the base node is selected on level 2, as indicated by the green checked circle on the left. The two grey checked fluid parcels labeled "a" and "d" are randomly selected. These two fluid parcels would then be swapped. This would change the pairing of parcels: \f$(a, b)\rightarrow(a,d)\f$ and \f$(c, d)\rightarrow (c,a)\f$. The second swapping process selects the right blue hash-filled circle at level 1 on the right. The two orange-hashed grandchild nodes are randomly selected, and the two subtrees below these nodes are swapped. In this case, those subtrees consist of fluid parcels \f$i\f$ and \f$j\f$ being swapped with parcels \f$m\f$ and \f$n\f$. Only the first case swaps can change particle adjacency (pairing) and thus influence parcels. Other swaps can affect the flow state but do not directly influence mixing. However, their role is crucial because they reduce the length scales (in the proximity sense) of parcel-to-parcel property variations and thus facilitate the eventual pairing of dissimilar parcels that induces molecular mixing. We can think of parcel mixing at level \f$N-3\f$ as micromixing and node/subtree mixing at other levels as macromixing that affects mixing length scales but not fine-scale fluid properties at the parcels.




## Eddy Selection

In the HiPS tree, the eddy rate \f$\lambda_i\f$ at level \f$i\f$ is determined as \f$1/\tau_i\f$ for each node, and collectively for all nodes at level \f$i\f$, it is calculated as:

$$
\lambda_i = \frac{2^i}{\tau_i}.
$$
The occurrence times of eddy events are sampled from an exponential distribution, representing a Poisson process with a mean rate \f$\Lambda\f$. This distribution is given by:

$$
p(\Delta t) = \Lambda e^{-\Lambda\Delta t}.
$$
To generate eddy occurrence times, we can sample from this distribution using:

$$
\Delta t = -\frac{\ln(P_r)}{\Lambda},
$$
where \f$P_r\f$ denotes a uniform random variate in the range \([0,1]\).

To determine the tree level for a sampled eddy event, we follow a process that captures the full range of scalar fluctuations from the inertial range of scales, \f$l \ge l^*\f$, denoted as \f$I\f$, to the viscous scale range, \f$l \le l^*\f$, denoted as \f$V\f$. Here, \f$l^*\f$ is the characteristic length scale delineating the transition between the inertial and viscous regimes, and \f$l\f$ is an arbitrary length scale. 
The total rate is \f$\Lambda = \Lambda_I + \Lambda_V\f$. A uniform random variate \f$P_r\f$ selects the region. If \f$P_r \le \Lambda_I/\Lambda\f$, region \f$I\f$ is chosen; otherwise, region \f$V\f$ is selected, and a specific level in the chosen region is then selected.

Within the inertial range, the probability of an eddy event at level \f$i\f$ is given by:

$$
p(i) = \frac{\lambda_i}{\Lambda_I},
$$

where \f$\Lambda_I = \sum_{i=0}^{i^*} \lambda_i\f$. The cumulative distribution function (CDF) is defined as:

$$
P(i) = \frac{1-(2A^{-2/3})^{i+1}}{1-(2A^{-2/3})^{i^*+1}},
$$

and level \f$i\f$ can be sampled as:

$$
i = \left\lceil \frac{\log_2( 1-P_r(1-(2A^{-2/3})^{i^*+1}) )}{1-\frac{2}{3}\log_2A} - 1 \right\rceil.
$$

Within the viscous range, the probability of an eddy event at level \f$i\f$ is given by:

$$
p(i) = \frac{\lambda_i}{\Lambda_V},
$$

where \f$\Lambda_V = \sum_{i=i^* +1}^{N-3}\lambda_i\f$. The CDF is expressed as:

$$
P(i) = \frac{2^{i+1}- 2^{i^* +1}}{2^{N-2}- 2^{i^* +1}},
$$

and level \f$i\f$ can be sampled using a uniform random variate \f$P_r\f$ as:

$$
i = \left\lceil\log_2\left(P_r(2^{N-2}-2^{i^* +1})+2^{i^* +1}\right) -1 \right\rceil.
$$


##  Schmidt Dependence

The Schmidt number (Sc) represents the ratio of momentum diffusivity to species diffusivity. In our context, we extend this definition to apply to arbitrary scalars, not limited to chemical species. The variable Sc formulation in the HiPS model varies depending on whether Sc is greater than 1 or less than 1. This discrepancy arises from differences in how mixing is handled in the HiPS model compared to diffusion-based simulation approaches that typically model transport processes in physical coordinates. We first present the variable Sc model regarding scalars aligned with tree levels in length and timescales before extending it to arbitrary scales.

The \f$l^*\f$ scale holds particular significance for a scalar with Sc equal to unity and can be directly compared to the Kolmogorov scale \f$\eta\f$, with the two scales being proportional. Additionally, \f$l^*_s\f$ represents the smallest length scale for a scalar with an arbitrary \f$Sc\f$, akin to the Batchelor scale \f$\eta_b\f$ (for \f$Sc\ge 1\f$) or the Obukhov-Corrsin scale \f$\eta_{oc}\f$ (for \f$Sc\le 1\f$), and is proportional to \f$\eta_b\f$ or \f$\eta_{oc}\f$, respectively. The HiPS Schmidt number is defined as:

$$
Sc=(l^*/l_s^*)^{p_s}. 
$$
\label{eq:ScDef}

For \f$Sc\ge 1\f$, the choice of \f$p_s=2\f$ is grounded in the analogy to a physical flow, where \f$\tau_\eta=\tau_{\eta_b}\f$ in the viscous regime. Here, \f$Sc=\nu/D\f$, \f$\nu=\eta^2/\tau_\eta\f$, and \f$D=\eta_b^2/\tau_{\eta_b}\f$.

For \f$Sc\le 1\f$, the selection of \f$p_s=3/4\f$ is similarly based on an analogy to a physical flow, where \f$\tau_{\eta_b}=\tau_\eta(\eta_b/\eta)^{2/3}\f$ in the inertial range (using \f$\tau_i = \tau_0A^{2i/3}\f$).

### Discrete

To establish the relationship between \f$Sc\f$ and the associated scales in the HiPS tree, we initially focus on \f$l^*_s\f$ values restricted to HiPS levels, followed in Sec.~\ref{s:arb_sc} by generalization to arbitrary \f$l^*_s\f$. In a general HiPS simulation, multiple scalars with different \f$Sc\f$ values can be considered.

The \f$Sc\f$ value is related to the tree levels \f$i^*\f$ and \f$i_s^*\f$ as follows:
$$
Sc = A^{p_s(i^*-i^*_s)} = A^{-p_s\Delta i}.
$$
\label{eq:Sci}

For \f$A=1/2\f$ and \f$p_s=4/3\f$ for \f$Sc\le 1\f$, we get \f$Sc=4^{2\Delta i/3}\f$, resulting in \f$Sc\approx\f$ 1, 0.4, 0.16, 0.062, 0.025, for \f$\Delta i=\f$ 0, -1, -2, -3, -4, respectively. For \f$Sc\ge 1\f$, with \f$p_s=2\f$, \f$Sc=4^{\Delta i}\f$, leading to \f$Sc=\f$ 1, 4, 16, 64, 256, for \f$\Delta i\f$ = 0, 1, 2, 3, 4, respectively.

### Arbitary


The formulation presented above identifies \f$Sc\f$ values associated with integer levels of the HiPS tree. Arbitrary \f$Sc\f$ corresponds to scalars with \f$l^*_s\f$ between two HiPS tree levels, and therefore, \f$i^*_s\f$ may not be an integer. For a scalar with a given \f$Sc\f$, \f$i^*_s\f$ is computed as \f$i^*_s = i^* - (\log Sc)/(p_s\log A)\f$, or

\f$
Sc<1: \phantom{xxx}    i^*_s = i^* - \frac{3\log Sc}{4\log A},
\f$
\f$
Sc>1: \phantom{xxx}    i^*_s = i^* - \frac{\log Sc}{2\log A}.
\f$
\label{eq:Sc<1}

Levels \f$i_-\f$ and \f$i_+\f$ are considered for the lower and upper levels bounding \f$i^*_s\f$. Eddy events that occur on levels at or above \f$i_+\f$ mix the scalar across the left and right subtrees emanating from the \f$i_+\f$-level eddy node. For a level \f$i_-\f$ eddy event, the scalar is mixed across the two subtrees of the \f$i_-\f$-level node with probability \f$p_-\f$, where

$$
p_- = \frac{i_+ - i^*_s}{i_+-i_-} = i_+-i^*_s.
$$
\label{e:p_-}

The second equality holds since \f$i_+-i_-\f$ is always unity. This probability is linear in index space and takes a value of 1 when \f$i^*_s=i_-\f$ and 0 when \f$i^*_s=i_+\f$. Using Eqs.~\ref{e:L}, \ref{e:tau}, and \ref{e: erate}, \f$p_-\f$ can be written as

$$
p_- = \frac{\log(l^*_s/l_+)}{\log(l_-/l_+)} = \frac{\log(\lambda^*_s/\lambda_+)}{\log(\lambda_-/\lambda_+)}.
$$

This form illustrates that while there is a linear interpolation of \f$p_-\f$ in index space, the corresponding interpolation between eddy lengths, times, or rates is logarithmic, consistent with the geometric progression of the scales with tree level.


