# Instruments

1. Cost shifters: cost shifters that vary across products
2. BLP instruments: functions of other product characteristics
3. Prices in other independent markets
4. Differentiation instruments: flexible approximations of an unknown function of the difference (quadratic and spline)
5. Chamberlain optimal instruments: 
6. Exogenous policy shocks: tax change
7. In the spirit of DiD

# Example

## Simple multinomial logit
$$\log(s_j)-\log(s_0)=x_j\beta+p_j\alpha+\xi_j$$
In simple logit multinomial model, we instrument the endogenous price with either the cost shifters or the BLP instruments.
In BLP1995, the first order basis functions associated with $z_{jk}$ 
1. $z_{jk}$: own product characteristics
2. $\sum_{r\neq j, r\in F} z_{rk}$: other product characteristics from the same firm
3. $\sum_{r\neq j, r\notin F} z_{rk}$: other product characteristics from rival firms

## Nested logit
$$\log(s_j)-\log(s_0)=x_j\beta+p_j\alpha+\sigma\log(\bar{s}_{j|g})+\xi_j$$
In nested logit model, we need more instruments to instrument the endogenous market share in the group $\log(\bar{s}_{j|g})$.
These variables might include the characteristics of other firms in the group.
1. $z_{jk}$: own product characteristics
2. $\sum_{r\neq j, r\in F} z_{rk}$: other product characteristics from the same firm in the group
3. $\sum_{r\neq j, r\notin F} z_{rk}$: other product characteristics from rival firms in the group