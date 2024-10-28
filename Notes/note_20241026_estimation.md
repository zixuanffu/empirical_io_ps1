# Simple multinomial logit

## Utility, probability and market share


$$ u_{ij}=\beta_1x_{j1}+\cdots+p_j\alpha+(\beta_0+\xi_j)+\epsilon_{ij}=\delta_j+\epsilon_{ij}$$
where the $(\beta_0+\xi_j)
$ is the product fixed effect. 

Assuming that $\epsilon_{ij}$ is T1EV and independent across all $j$ (we do not assume independence across all $i$!),**we will relax this assumption and capture the dependce explicitly later). We have
$$ P_{ij}=\frac{e^{u_{ij}}}{e^{u_{i0}}+\sum_{k=1}^{J}e^{u_{ik}}} $$
Since now, we assume all individuals have identical coefficients $\beta,\alpha$.The market share of product $j$ is the choice probability of product.

## Estimation
we normalize the utility of the outside option to zero, so the market share $s_j$ can be written as 
$$ s_j=\frac{e^{u_{ij}}}{1+\sum_{k=1}^{J}e^{u_{ik}}} $$
Then the log-odds ratio
$$\log(\frac{s_j}{s_0})=x_j\beta+p_j\alpha+\xi_j$$
that is
$$ \log s_j-\log s_0=\beta_0+\beta_1x_{j1}+\cdots+p_j\alpha+\xi_j$$
while the $\xi_j$ is the *error term* in a linear model. In order to use the OLS regression, we need to make the error term independent across all $j$. 

However, this is very unlikely because product fixed effect $\xi_j$ can correlat with price.
> In panel data, the individual FE is constant across time and we can take difference to eliminate it. Here, if there are multipe markets markets where the set of products is sold we may be able to take difference (?)

# Multinomial logit with nests

Previously we assume that $\epsilon_{ij}$ is T1EV and independent across all $j$, now we capture the dependence across all products by **nesting**. That is to say, the utility shocks are correlated for those within the same nest.

## Utility, probability and market share
$$u_{ij}=x_j\beta+p_j\alpha+\xi_j+\zeta_{ig}+(1-\sigma)\epsilon_{ij}=\delta_j+\zeta_{ig}+(1-\sigma)\epsilon_{ij}$$


## Estimation

$$\log(s_j)-\log(s_0)=x_j\beta+p_j\alpha+\sigma\log(\bar{s}_{j|g})+\xi_j$$
where $\bar{s}_{j|g}$ is the average market share of products in the same nest.

# Multinomial logit with interacted taste coefficients

## Utility
$$u_{ij}=x_j\beta z_i+p_j\alpha z_i+(\beta_0+\xi_j)+\epsilon_{ij}$$
where $z_i$ is the individual characteristics. such as income, age. 
Then 
$$u_{ij}=\delta_j+x_j\beta (z_i-\bar{z})+p_j\alpha(z_i-\bar{z})+\epsilon_{ij}$$
where $\bar{z}$ is the average of $z_i$.
Terms after $\delta_j$ are to capture the dependence of utility shocks across all products.

# Multinoimal logit with random coefficients

This is a generalization of the previous model. Instead of having taste coefficients specific to income/age, we allow it to be even more flexible, specific to each individual!
The purpose is also to capture the dependence of utility shocks across all products, just to reemphasize.

## Utility
$$u_{ij}=x_j\beta_i+p_j\alpha_i+(\beta_0+\xi_j)+\epsilon_{ij}= \delta_j+\sum_k \sigma_k x_{jk}v_{ik}+\epsilon_{ij}$$

$$ u_{ij}=\delta_j+\mu_{ij}$$
where $\mu_{ij}$ captures the dependence.

where $\beta_i=\bar{\beta}+\sigma v_i$ and $\alpha_i=\bar{\alpha}+\sigma w_i$ with $v_i$ follows $N(0,\Sigma_v)$. 

# Side words

##  Log*
I want to clarify some terminology associated with *logit*.

### logistic distribution
The standard $(\mu=0,s=1)$ is 
$$ F(x)=\frac{1}{1+e^{-(x-\mu)/s}} $$

### logistic function
The (standard) logistic function is defined as the cdf of the logistic distribution.
$$\sigma(x)=\frac{1}{1+e^{-x}}$$

### logit function (log-odds)

$$\text{logit}(p)=\log(\frac{p}{1-p})$$

> The logit function is the quantile function of the standard logistic distribution.
> $$ F_X(x)=p \Rightarrow Q(p)=F_X^{-1}(p)$$

> Probit follows the same logic, with the probit function as the quantile function of the standard normal distribution. The *probistic* function is the cdf of the standard normal distribution.