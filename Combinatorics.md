# Combinatorics



## Permutation
$P^m_n = A^m_n = n(n-1)(n-2)\dots(n-m+1) = \dfrac{n!}{(n-m)!},\ m\leq n$

## Binomial
$\dbinom{n}{m} = \dfrac{P^m_n}{P_m} = \dfrac{n!}{m!(n-m)!},\ \dbinom{n}{0} = 1$

$\dbinom{n}{m} = \dfrac{n(n-1)(n-2)(n-3)\dots(n-m+1)}{m!},\ m\leq n$

## 有重复的排列

$从N个元素中取出M个元素(有放回)的排列方案数为N^{M}$

## 有重复的组合

$从N个元素中取出M个元素(有放回)的组合方案数为\dbinom{n + m - 1}{m}$

### 隔板法

$m个球,\ 放到n个盒子里面(允许空盒子)的方案数$

$\ x_i表示第i个盒子中球的个数.\ \sum\limits^{n}_{i = 1}{x_i} = m, x_i \geq 0. 相当于求左式中方程解的个数$

$这个问题等价于m个球中间插入n - 1个板, 板与板之间允许为空$

$等价于m + n - 1个位置中选择n - 1个位置放板子, 方案数为\dbinom{m + n - 1}{n - 1} = \dbinom{m + n - 1}{m}$



$如果问题改为x_i > 0$

$可以看作x_i\prime = x_i - 1, x_i\prime \geq 0, \sum\limits^{n}_{i = 1}{x_i\prime} = m - n$

$答案为\dbinom{m - n + n - 1}{n - 1} = \dbinom{m - 1}{n - 1}$

$也可以直接看有m - 1个位置, 有n - 1种选法, 那么答案就是\dbinom{m - 1}{n - 1}$



### 有相同类型元素的排列数

#### 比如3个红球, 2个绿球, 2个🏀的不同排列数

$$
\binom{7}{3} \binom{4}{2} \binom{2}{2} = \dfrac{7!}{3!\ 2!\ 2!}
$$

$推广: 有n个位置, 要填a_i 个i, 1 \leq i \leq m$
$$
答案: \dbinom{n}{a_1}\dbinom{n-a_1}{a_2}\dbinom{n-a_1-a_2}{a_3}\dots\dbinom{n - \sum\limits^{k - 1}_{i = 1}{a_i}}{a_k}\dots\dbinom{a_m}{a_m} = \dfrac{n!}{\prod\limits^{m}_{i = 1}{(a_i!)}}
$$




## 二项式定理

$$
(x + y)^n = \sum\limits^{n}_{k = 0}{\dbinom{n}{k}x^{n-k}y^k} = \sum\limits^{n}_{k = 0}{\dbinom{n}{k}x^{k}y^{n - k}}
$$

组合意义上:

$有n个(x + y)相乘, 每次都从(x + y)中选择一个x或者选择一个y, 那么多项式某一项的形式形如x^{k}y^{n - k}(n个式子中选了k个x, n - k个y)$
$多项式的系数为\dbinom{n}{k}表示为从n种选择中选了k个x的方案数$

特别的:
$$
(x + 1) ^ {n} = \sum\limits^{n}_{i = 0}{\binom{n}{i}x^i}
$$


## 一些恒等式

### 基本恒等式

$$
k\dbinom{n}{k} = n\dbinom{n - 1}{k - 1}
$$

$$
\dbinom{n}{k}\dbinom{m}{k} = \dbinom{m}{n}\dbinom{m - n}{m - k} (要求:m - k < m - n)
$$

$$
\sum\limits^{n}_{i = 0}{\binom{n}{i}} = 2^n
$$

特别的:
$$
\sum\limits^{n}_{\substack{i = 0\\ i \equiv 0 \pmod{2}} }{\binom{n}{i}} = \sum\limits^{n}_{\substack{i = 1\\ i \equiv 1 \pmod{2}} }{\binom{n}{i}} =2^{n-1}
$$




### $(1 + x)^n$二项式展开

$$
\sum\limits^{n}_{k = 0}{(-1) ^ {k}\binom{n}{k}} = 0
$$

$$
\binom{n}{k} + \binom{n}{k + 1} = \binom{n + 1}{k + 1}
$$

### 杨辉三角推广

固定选的元素个数求和, 可以优化速度
$$
\sum\limits^{m}_{k = 0}{\binom{n + k}{k}} = \sum\limits^{m}_{k = 0}{\binom{n + k}{n}} = \binom{n + m + 1}{m + 1} = \binom{n + m + 1}{n}
$$

## 求幂和

$$
\sum\limits^{n}_{i = 1}{i} = \sum\limits^{n}_{i = 1}{\binom{i}{1}} = \binom{n + 1}{2} = \frac{n(n + 1)}{2}
$$

------

高次写成下降幂: 
$$
i^2 = i(i-1) + i = 2\binom{i}{2} + \binom{i}{1}
$$

$$
\sum\limits^{n}_{i = 1}{i^2} = \sum\limits^{n}_{i = 1}{i(i-1)+i} = 2\binom{n + 1}{3} + \binom{n + 1}{2} = \frac{n(n + 1)(2n + 1)}{6}
$$


------


$$
i^3 = (i+1)i(i-1) + i = 6\binom{i+1}{3} + \binom{i}{1}
$$

$$
\sum\limits^{n}_{i = 1}{i^3} = \sum\limits^{n}_{i = 1}{(i + 1)i(i-1)+i} = 6\binom{n + 2}{4} + \binom{n + 1}{2} = \frac{n^{2}(n+1)^{2}}{4} = \binom{n + 1}{2} ^ 2
$$

## 整数划分

把一个正整数$N$写作多个大于等于1且小于等于其本身的整数的和, 则其中各个加数所构成的集合为$N$的一个划分

看作一个多重背包, 背包的体积为$N$, 有$N$种物品, 体积分别为$1\dots N$, 恰好把这个多重背包填满的方案数

增长的没有想象的那么快

## 第二类斯特林数(集合拆分)

$把N个元素划分到K个非空集合中, 有多少种做法$

$S(n, m) 为方案数$
$$
S(n, 0) = 0, \ (k = 0)\\
S(n, 1) = 1, \ (k = 1)\\
S(n, n) = 1\\
S(n, k) = 0, \ (k > n)\\
S(n, k) = kS(n - 1, k) + S(n - 1, k - 1), \ (1 < k < n)
$$
前四个是初值和边界情况

最后一个是枚举当前第$n$个元素是否放到新的集合第$k$个集合中

$S(n - 1, k - 1)是把n单独划分集合的情况, 即把n - 1个元素划分到前k - 1个集合中, 第n个元素划分到第k个集合中的情况$

$kS(n - 1, k) 是把第n个元素放到前面k个中的集合中的某一个的情况, 有k种不同放法, 所以前面系数要乘k$

##  卡特兰数

一个数列: $\{C_n\} = 1, 1, 2, 5, 14, 42, 132, \dots$

1. 有一个大小为$N * N$的方格图左下角$(0, 0)$为右上角为$(n, n)$，从左下角开始每次都只能向右或者向上走一单位，不走到对角线$y=x$上方（但可以触碰)的情况下到达右上角有多少可能的路径
2. 在圆上选择 $2 * N$ 个点，将这些点成对连接起来使得所得到的 $N$ 条线段不相交的方法数
3. 对角线不相交的情况下，将一个凸多边形区域分成三角形区域的方法数
4. 一个栈（无穷大）的进栈序列为$1, 2, 3, \dots, n$有多少个不同的出栈序列？
5. $N$个结点可构造多少个不同的二叉树？

满足如下公式的
$$
C_0 = C_1 = 1\\
C_i = \sum\limits^{n-1}_{i = 0}{C_{i}C_{n-1-i}}, \ i \geq2
$$
其他公式
$$
C_n= \frac{(4n-2)C_{n-1}}{n + 1}\\
C_n = \binom{2n}{n} - \binom{2n}{n-1}
$$
