# Supply side 

When firms change ownership structure while keeping the marginal cost constant, firms will change price because now that they face a new profit maximisation problem.

In order to calculate the new price, we need to know the existing marginal cost. 

## Without modeling the marginal cost 
This section is devoted to backing out the marginal cost from the observed price and market share.

The profit of a multiproduct firm is 
$$ \Pi_f = \sum_{j\in \mathcal{J}_f} (p_j -c_j )s_jM$$
since $s_jM$ represents the sales of the product $j$.

The optimal price is given by taking the first order derivative wrt p.

$$ s_j+\sum_{k\in \mathcal{J}_f} \frac{\partial s_k}{\partial p_j} (p_k-c_k) = 0 \quad \forall j\in \mathcal{J}_f$$

The first term is the direct impact of changing $p_j$. The second term is the impact on the profits of other products of the firm.

For firm $f$, the pricing equation written in the matrix form gives,
$$ S^f + \Omega^f(P^f-C^f) = 0$$
where $\Omega^f$ is the matrix of cross-price elasticities. For example, if firm $f$ has two products, the matrix is
$$ \Omega^f =\begin{bmatrix} 
\frac{\partial s_1}{\partial p_1} & \frac{\partial s_1}{\partial p_2} &0 \\\
\frac{\partial s_2}{\partial p_1} &\frac{\partial s_2}{\partial p_2} & 0 \\0
& 0 &0\\
\end{bmatrix}$$

Therefore, aggregate the pricing equation for all firms, we have
$$ S + \Omega \odot O (P-C) = 0 $$
where $O$ is the ownership matrix.

The ownership matrix is computed by first create a dummy matrix $D$ of dimension $J\times F$ (the number of products times the number of firms). Then $JJ^T$ is the ownership matrix $O$.

## Modelling the marginal cost with cost shifters.

$$ c_j = w_j \gamma + \omega_j$$

1. if marginal cost depends on quantity
2. predict marginal cost of new product
3. gain precision on the estimation of demand parameters. (Instrument for the price?)