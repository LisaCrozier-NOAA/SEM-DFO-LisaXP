# Here is exactly where the math breaks down and why `lavaan` is giving you **$df = 3$**.
# 
# It comes down to a fundamental rule of structural equation modeling: **Degrees of freedom have absolutely nothing to do with your sample size ($N = 24$ or $N = 23$).** Sample size affects your power and standard errors, but $df$ is calculated strictly using the number of variables, the number of estimated parameters, and—critically—**whether your model includes a mean structure.**
#   
#   Because you have a model with regressions, intercepts, and missing data, `lavaan` automatically switches on a **mean structure** to compute Full Information Maximum Likelihood (FIML). This changes the information math completely.
# 
# ---
#   
#   ### 1. The True Total Information Available
#   
#   When a mean structure is included, R doesn't just count the unique variances and covariances $\frac{p(p+1)}{2}$. It **also** includes the mean/intercept of every single variable ($p$).
# 
# The corrected formula for total pieces of information is:
# 
# 
# $$\text{Total Information} = \frac{p(p+1)}{2} + p$$
# 
# You mentioned your model has **7 variables** (nodes). Let's plug 7 into the formula:
#   
#   
#   $$\text{Total Information} = \frac{7(7+1)}{2} + 7 = 28 + 7 = 35$$
#   
#   So, your data starts with **35** total pieces of available information.
# 
# ---
#   
#   ### 2. The Free Parameters Count
#   
#   Your output explicitly states: `Number of model parameters = 12`.
# 
# However, because this is an SEM path model, `lavaan` is also automatically estimating the variances and means of your **5 exogenous predictors** behind the scenes to handle the missing data patterns.
# 
# Let's break down the *actual* total number of parameters being estimated across your 7 variables:
# 
# * **8** Regression Coefficients (Explicitly in your output)
# * **2** Intercepts for your 2 dependent variables (Explicitly in your output)
# * **2** Residual Variances for your 2 dependent variables (Explicitly in your output)
# * **5** Means/Intercepts for your 5 exogenous predictor variables (Implicitly estimated by FIML)
# * **15** Variances and Covariances among your 5 exogenous predictor variables ($\frac{5 \times 6}{2} = 15$, implicitly estimated)
# 
# $$\text{Total Estimated Parameters} = 8 + 2 + 2 + 5 + 15 = 32$$
# 
# ---
# 
# ### 3. The Final $df$ Calculation
# 
# Now we can subtract our total parameters from our total available information pool:
# 
# $$df = \text{Total Information} - \text{Total Estimated Parameters}$$
# 
# $$df = 35 - 32 = 3$$
# 
# This is exactly how `lavaan` arrives at **3 degrees of freedom**.
# 
# Your model leaves exactly 3 possible paths or relationships among your variables completely unestimated (fixed to zero), allowing the software to test your model fit with a Chi-Square statistic!
# 
# ---
# 
# For an excelent visual breakdown of how `lavaan` computes and reports its fit statistics under FIML, you can review this [lavaan Missing Data Tutorial Video](https://www.youtube.com/watch?v=ofzZA9J5-eI) which covers how degrees of freedom are preserved and evaluated when handling missing datasets.