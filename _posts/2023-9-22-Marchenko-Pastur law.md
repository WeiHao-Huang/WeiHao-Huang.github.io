>  本文旨在以简短的篇幅阐述 Marchenko-Pastur Law 及其相关的应用。

<!-- TODO:Marchenko-Pastur Law 要不要打法语音标？-->

在许多实际问题（例如：核物理以及量子力学）中需要得到大维随机矩阵（Large Dimensional Random Matrices，LDRM）的谱。

然而，当矩阵维度较高的时候，经验谱分布（Empirical Spectral Distribution，ESD）的形式非常复杂，因此，人们将注意力转向极限谱分布（Limiting Spectral Distribution，LSD）。其中一个著名的结果是 Marcenko-Pastur law，它给出了大维样本协方差矩阵（large dimensional sample covariance matrices）的极限谱分布。

本文讨论如下几个主题：极限谱分布的存在性、LSD 的显式表达式、ESD 的最大最小特征值、 ESD 到 LSD 的收敛速度。



<h1>Preliminary</h1>

<h2>Sample covariance matrices</h2>

<!--全文仅仅考虑实矩阵要不要说？-->

设 $\boldsymbol{X}_{1}, \boldsymbol{X}_{2}, \dots, \boldsymbol{X}_{n}$ 是独立同分布的 $m$ 维随机列向量，其中，$\boldsymbol{X}_{k} = (x_{1k}, \dots, x_{mk})^{{\rm T}}$ 满足均值为 $\boldsymbol{0}$，协方差矩阵为 $\sigma^{2} I_{m}$。如此可以定义 $X = ( \boldsymbol{X}_{1}, \boldsymbol{X}_{2}, \dots, \boldsymbol{X}_{n} ) \in \mathbb{R}^{m\times n}$ 是 $m \times n$ 维的随机矩阵。

现在可以来定义研究对象：

**定义**：样本协方差矩阵 $S$
$$
\begin{equation}
	S = \frac{1}{n-1} \sum_{k=1}^{n} ( \boldsymbol{X}_{k}-\overline{\boldsymbol{X}} )( \boldsymbol{X}_{k}-\overline{\boldsymbol{X}} )^{\rm T}
\end{equation}
$$
其中，$\overline{\boldsymbol{X}} = \frac{1}{n} \sum_{k=1}^{n} \boldsymbol{X}_{k}$。

但在大维随机矩阵中，样本协方差矩阵通常定义为 $S = \frac{1}{n} \sum_{k=1}^{n} \boldsymbol{X}_{k} \boldsymbol{X}_{k}^{\rm T} = \frac{1}{n} X X^{\rm T} \in \mathbb{R}^{m \times m}$。因为不减去样本均值向量并不会影响 $S$ 的极限谱分布，这可以通过 rank inequality 证明（Bai 于 1999 年给出了此证明），本文接下来的内容均采用这种定义。



<h2>ESD</h2>

为了得到极限谱分布，需要先定义经验谱分布。

为此，先设 $S$ 的特征值为 $\lambda_{1}, \lambda_{2}, \dots, \lambda_{m}$；于是可以定义如下的随机测度：
$$
\begin{equation}
	\mu_{m} = \frac{1}{m} \sum_{i = 1}^{m} \delta_{\lambda_{i}}
\end{equation}
$$
该随机测度计算 $\mathbb{R}$ 的子集 $A$ 中特征值的数量 $\mu_{m}(A)$。

**定义**：经验谱分布 $F^{S}$
$$
\begin{equation}
	F^{S}(x) = \mu_{m}((-\infty, x])
\end{equation}
$$
在此基础上，Marcenko 和 Pastur 于 1967 年首次得出了 $S$ 的极限谱分布，我们现在来了解此分布的形式。



<h1>Marcenko-Pastur law</h1>



这一部分旨在介绍 Marcenko 和 Pastur 得出的关于 $S$ 的极限谱分布，现在这个分布被称为 Marcenko-Pastur law。

关于 $S$ 的 LSD，有四个问题值得讨论：LSD 的存在性、 LSD 的显式表达式、极限特征值、以及 ESD 收敛于 LSD 的速度。本节剩余的内容就来讨论这四个主题。



<h2>Existence of LSD</h2>

结合 Moment Convergence Theorem (MCT) 和 Carleman’s condition 可以证明经验谱分布 $F^{S}$ 收敛于极限谱分布 $F$（更一般的说，对任意 Hermitian 矩阵都有极限谱分布存在的结论）。在极限谱分布存在的情况下，只有一部分矩阵可以被确定出分布的显式表达式，幸好，大样本协方差矩阵 $S$ 就存在分布的显式表达式。



<h2>Explicit expression</h2>

**定理**

假设 $m,n \to \infty$，且 $\lim_{m,n \to \infty} \frac{m}{n} = \rho \in (0, +\infty)$，那么随机测度 $\mu_{m}$ 收敛（weak* topology convergence）到一个确定的测度 $\mu$（$\mu$ 是 $F$ 的概率密度）：
$$
\begin{equation}
	\mu = \left\{
	\begin{aligned}
		& (1-\frac{1}{\rho}) \delta_{0} + \nu, & \rho >1\\
		& \nu, & 0 < \rho \leq 1
	\end{aligned}
	\right.
\end{equation}
$$
$\nu$ 是一个测度，其 Radon-Nikodym 导数满足：
$$
\begin{equation}
	\frac{{\rm d} \nu}{{\rm d}x} = \frac{1}{2\pi \sigma^{2}} \frac{\sqrt{(\lambda_{+} - x)(x - \lambda_{-})}}{\rho x} 
\end{equation}
$$
其中，
$$
\begin{equation}
	\lambda_{-} = \sigma^{2} (1 - \sqrt{\rho})^{2}, \quad \lambda_{+} = \sigma^{2} (1 + \sqrt{\rho})^{2}
\end{equation}
$$
$\sigma$ 是尺度参数，$\rho$ 是比例系数，接下来我们讨论比例系数 $\rho$ 对极限谱分布的影响。



当 $\rho > 1$ 时，极限谱分布为
$$
\begin{equation}
	(1 - \frac{1}{\rho}) \delta_{0} + \nu
\end{equation}
$$
此时，该概率由纯不连续的部分 $(1 - \frac{1}{\rho}) \delta_{0}$ 和连续的部分 $\nu$ 组成。

若 $\rho \leq 1$，那么有
$$
\begin{equation}
	\frac{\mathrm{d} \mu}{\mathrm{d} x} = \frac{1}{2\pi \sigma^{2}} \frac{\sqrt{(\lambda_{+} - x)(x - \lambda_{-})}}{\rho x}
\end{equation}
$$


<h2> Limits of extreme eigenvalues</h2>

即便有了 LSD 的显示表达式，但是我们对大样本协方差矩阵的最大、最小特征值的分布依然知之甚少，于是在 Geman, Yin, Bai, Krishnaiah, Silverstein 等人先后努力下，最大、最小特征值的强收敛形式被建立起来。

具体来说，在随机矩阵 $S$ 每一个元素都存在有限四阶矩的条件下，我们可以得到：
$$
\begin{equation}
	\left\{
	\begin{aligned}
		& \lim_{n \to \infty} {\rm min}\; \lambda(S) = \sigma^{2} (1 - \sqrt{\rho})^{2}, \; a.s. \\
		& \lim_{n \to \infty} {\rm max}\; \lambda(S) = \sigma^{2} (1 + \sqrt{\rho})^{2}, \; a.s. 
	\end{aligned}
	\right.
\end{equation}
$$
该结果证明了在极限情况下，$S$ 的最大最小特征值消除了随机性，具有确定的值。

有趣的是，通过统计物理中的 replica method 也可以导出以上的结论，具体内容有待补充。

<!--COMMENT: SPECTRAL ANALYSIS OF RANDOM MATRICES USING THE REPLICA METHOD Wishart note. To be added.-->



<h2>Convergence rate of the ESD</h2>

Bai 等人于 1993 年证明了，如果以下条件被满足：

- $\mathbb{E}[x_{jk}] = 0$，对所有的 $j, k, n$；
- $\mathbb{E}[x_{jk}^{2}] = 1$，对所有的 $j, k, n$；
- $ \sup_{n} \sup_{j, k} \mathbb{E}[|x_{j, k}|^{4}] I_{(|x_{j, k}| \geq M)} \to 0; M \to \infty $

那么下述结论成立：



**定理**

当 $0 < \theta < \Theta <1$ 或 $1 < \theta < \Theta < \infty$ 时，有：
$$
\begin{equation}
	\mathop{\sup}\limits_{\rho \in (\theta, \Theta)} \Vert \mathbb{E}[F^{S}] - F_{\lambda} \Vert = O(n^{-\frac{1}{4}})
\end{equation}
$$
**定理**

当 $0 < \varepsilon < 1$ 时，有：
$$
\begin{equation}
	\mathop{\sup}\limits_{\rho \in (1-\varepsilon, 1+\varepsilon)} \Vert \mathbb{E}[F^{S}] - F_{\lambda} \Vert = O(n^{-\frac{5}{48}})
\end{equation}
$$


Bai，Miao 以及 Tsay 等人于 1996 年推广了上述结论：

**定理**
$$
\begin{equation}
	\mathop{\sup}\limits_{\rho \in (\theta, \Theta)} \Vert F^{S} - F_{\lambda} \Vert = O_{p}(n^{-\frac{1}{4}})
\end{equation}
$$
**定理**
$$
\begin{equation}
	\mathop{\sup}\limits_{\rho \in (1-\varepsilon, 1+\varepsilon)} \Vert F^{S} - F_{\lambda} \Vert = O_{p}(n^{-\frac{5}{48}})
\end{equation}
$$



<h1>Application</h1>

这一部分旨在讨论 Marcenko-Pastur law 的应用： To be added.



<!--Expand: To relax the independence of the entries of Xk：LARGE SAMPLE COVARIANCE MATRICES WITHOUT INDEPENDENCE STRUCTURES IN COLUMNS-->

<!--google search:the application of Marchenko-Pastur law-->
