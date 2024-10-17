The codes are designed for the simulation study.

# Setup

We considered a mixed linear model $$y=Xb+e=X(\beta+u)+e\in\\mathbb{R}^5$$ where $X=((1,1,1,1,1)^\intercal, (0,1,2,3,4)^\intercal$ and $b$ represent the random intercept and slope, which satisfy $u\sim N(0, \Sigma)$ and  $e\sim N(0, \sigma_e^2I)$ where

$$
\Sigma=
\begin{pmatrix}
\sigma_1^2 & \sigma_{12} \\
\sigma_{12} & \sigma_2^2
\end{pmatrix}\ .
$$

Let the parameters of interest be $\theta=(\sigma_1, \sigma_2, \sigma_{12}, \sigma_e, \beta_1, \beta_2)$. Eight response patterns are generated. The parent relations are randomly generated and remain fixed throughout the entire simulation. The mixture coefficients are all chosen as the first type. The propensity odds $$O^r(l^r)$$ are designed as functions of observed variables. The identifying assumption can be represented by the following graph.

<p align="center">
<img src="https://github.com/user-attachments/assets/d2e31d26-bfb4-4491-9305-7fc411b14022" />
</p>


# Result of an example


We assume that the identifying assumptions are correctly specified, meaning the pattern graph structure and choice of mixture coefficients are known. Initially, we analyzed the ideal scenario where the full simulated dataset (Full dataset) is observed. Subsequently, we excluded all missing data and analyzed only complete cases. Next, we employed inverse propensity weighting methods, starting with the true inverse propensity weights. These three approaches serve as benchmark comparisons for IPW methods using estimated weights. We implemented the aforementioned three approaches: ``Logistic`` denotes that the estimated propensity odds are transformations of estimated probabilities from standard logistic regressions. ``Separate`` and ``Sequential`` refer to two distinct balancing procedures, respectively. For a fair comparison, we applied the same basis functions and penalization strategy. The biases and mean squared errors are detailed in the following tables.

<p align="center">
<img src="https://github.com/user-attachments/assets/11947ac4-4367-46a2-b504-cef16219b666" />
</p>
