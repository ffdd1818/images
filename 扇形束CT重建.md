# 等角度扇形束（Equal-Angle Fan Beam）CT


 > !24 扇形束的数据足量条件：每一条与物体相交的直线都必须包含（至少）一个扇形束的焦点位置!
 > 
>  !24 任何直线，只要它与扇形的焦点轨迹有两个交点，那么沿这条线的线积分就被测到了两次!

![图片3.png](https://docimg4.docs.qq.com/image/AgAABenPl4FL_Z1UP-RPvJsUyI5G2gGB.png?w=1140&h=790)

## !!#ff0000 短扫描：!!当扇形束探测器旋转$2\pi$, 每条投影射线也都被测到了两次 ( 图 3.8) ，冗余的数据可由下面这个表达式给出：

![短扫描.png](https://docimg9.docs.qq.com/image/AgAABenPl4HqWPCkovxFcKRYlq00-PQg.png?w=775&h=604)

## 典型案例

**例题1**   下面图示的扇形束数据采集方案能为图像重建测得足够的数据吗? 这里，扇形束的焦点轨迹由三段不相连的圆弧组成 ( 图 3.13) 。
![图片4.png](https://docimg8.docs.qq.com/image/AgAABenPl4EMiE3Tk-BPv56Dk-CeX9V6.png?w=528&h=372)
<font color=red>解：</font>因为如果你过那个圆形物体随意画一条直线，这条线总会与扇形束的焦点轨迹至少有一点相交。

## 1) 重建算法原理
设射线源$S_o$到原点$O$（即扇面旋转中心）的距离为$D$，中心线左右扇面张开的角度为$\gamma_m$。**扇形位置**由该中心射线与y轴交角$\beta$(**旋转角度**)确定。因此在$x-y$坐标系中，射线的绝对位置由$(\beta, \gamma)$唯一确定。射线投影记为$p^f(\gamma, \beta)$，重建物体图像$f(x, y)$，用极坐标表示$f(r, \phi)$。如果把$S_oE$看作一条平行束投影系统中的射线，射线也可用$(t, \theta)$确定。

![图片1.png](https://docimg6.docs.qq.com/image/AgAABenPl4GhKT0YcfhGj7c3EJO_5O-8.png?w=801&h=543)

由图中的几何关系可知，
$$
\begin{array}{l}
\theta  = \beta  + \gamma \\
t = D\sin \gamma 
\end{array}
$$

已知平行束重建算法为（**卷积形式**）：
$$
f(x,y) = \int_0^\pi  {\int_{ - {t_m}}^{{t_m}} {p(t,\theta )h(x\cos \theta  + y\sin \theta  - t)dt} d\theta }
$$

其中$h(\cdot)$是滤波器。为了方便推导，换为极坐标形式，令$x=r \cos \phi, y=r\sin \phi$，则$x\cos \theta+y\sin \theta=r\cos(\theta-\phi)$，于是：
$$
f(r,\phi)=\dfrac{1}{2}\int_0^{2\pi}{\int_{ - {t_m}}^{{t_m}} {p(t,\theta )h[r\cos(\theta-\phi)  - t]dt} d\theta}
$$
利用$(t, \theta)\Leftrightarrow(\gamma, \beta)$之间的关系，将平行束的变量替换成扇形束的变量，替换后的微分关系可以利用雅克比行列式表示：
$$
\left| J \right| = \left| {\begin{array}{*{20}{c}}
{\frac{{\partial t}}{{\partial \gamma }}}&{\frac{{\partial t}}{{\partial \beta }}}\\
{\frac{{\partial \theta }}{{\partial \gamma }}}&{\frac{{\partial \theta }}{{\partial \beta }}}
\end{array}} \right| = \left| {\begin{array}{*{20}{c}}
{D\cos \gamma }&0\\
1&1
\end{array}} \right| = D\cos \gamma
$$
即$dtd\theta=D\cos \gamma d\gamma d\beta$，于是可得：
$$
f(r,\phi)=\dfrac{1}{2}\int_{-\gamma}^{2\pi -\gamma}{\int_{ - {arcsin(t_m/D)}}^{{arcsin(t_m/D)}} {p(D\sin\gamma,\beta+\gamma )\cdot h[r\cos(\beta+\gamma-\phi)  - D\sin\gamma]\cdot D\cos\gamma d\gamma} d\beta}
$$
考虑到$-\gamma \sim 2\pi-\gamma $已经覆盖了物体的整个$360^\circ$，所以等效于从0到$2\pi$进行积分，且$p(D\sin\gamma,\beta+\gamma )=p^f(\gamma,\beta)$，因此可得扇形束图像重建算法为：

$$
% \boxed{
% f(r,\phi)=\dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - \gamma_m}^{\gamma_m} {p^f(\gamma,\beta)\cdot
% h[r\cos(\beta+\gamma-\varphi)  - D\sin\gamma]\cdot D\cos\gamma d\gamma} d\beta}
%}
 \bbox[#f9de9f,6pt]{
 f(r,\phi)=\dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - \gamma_m}^{\gamma_m} {p^f(\gamma,\beta)\cdot
 h[r\cos(\beta+\gamma-\phi)  - D\sin\gamma]\cdot D\cos\gamma d\gamma} d\beta}
 }  
$$
**但上述公式对$\gamma$的积分不是!!#e60000 卷积!!的形式。**
为了配合卷积形式，根据下图对相关变量重新命名：$M$为重建点，重新定义射线夹角为$\gamma^\prime$，射线源到重建点的长度为$L$。

![图片2.png](https://docimg7.docs.qq.com/image/AgAABenPl4Gr264ow5BMiLRFlqTbG4RB.png?w=793&h=597)

$$
\begin{align}
   & r\cos(\beta+\gamma-\phi)-D\sin\gamma  \\
= & {\color{red}r\cos(\beta-\phi)}\cos \gamma-r\sin(\beta-\phi)\sin\gamma-D\sin\gamma \\
= & {\color{red}L\sin\gamma^\prime}\cos \gamma-{\color{blue}[r\sin(\beta-\phi)+D]}\sin\gamma  \\
= & {L\sin\gamma^\prime}\cos \gamma-{\color{blue} L\cos\gamma^\prime}\sin\gamma  \\
= & {L\sin(\gamma^\prime - \gamma)}
\end{align}
$$
于是扇形束图像重建算法可改写为：
$$
 f(r,\phi)=\dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - \gamma_m}^{\gamma_m} {p^f(\gamma,\beta)\cdot
 h[L\sin(\gamma^\prime - \gamma)]\cdot D\cos\gamma d\gamma} d\beta}
$$
> 斜坡滤波器卷积核的性质：$h[g(t)]=\dfrac{t^2}{g(t)^2}h(t)$，该性质只有在频域积分区间为无限时才成立。若对斜坡滤波器进行加窗，变成有限带宽的函数，则该等式不成立，会在重建的图像中造成不均匀的分辨率。

由该性质可得$h[L\sin( \gamma)]=\dfrac{\gamma^2}{[L\sin(\gamma)]^2}h(\gamma) $，进一步可得 $h[L\sin(\gamma^\prime - \gamma)]=\dfrac{(\gamma^\prime-\gamma)^2}{(L\sin(\gamma^\prime-\gamma))^2}h(\gamma^\prime-\gamma)$，于是扇形束FBP重建算法的卷积形式为：
$$
 \bbox[#f9de9f,6pt]{
\begin{align}
  f(r,\phi) &=  \dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - \gamma_m}^{\gamma_m} {p^f(\gamma,\beta)\cdot
 h[L\sin(\gamma^\prime - \gamma)]\cdot D\cos\gamma d\gamma} d\beta} \\
             &= \dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - \gamma_m}^{\gamma_m} {p^f(\gamma,\beta)\cdot
 \dfrac{(\gamma^\prime-\gamma)^2}{(L\sin(\gamma^\prime-\gamma))^2}h(\gamma^\prime-\gamma)\cdot D\cos\gamma d\gamma} d\beta}   \\
            &=\color{red} \int_{0}^{2\pi}{\dfrac{1}{L^2}[p^f(\gamma,\beta)D\cos\gamma]*\dfrac{\gamma^2}{2\sin^2\gamma}h(\gamma)d\beta}
\end{align}
}
$$

其中：
$$
\begin{align}
  L&=\sqrt{(D+r\sin(\beta-\phi))^2+(r\cos(\beta-\phi))^2} =\sqrt{D^2+r^2+2Dr\sin(\beta-\phi)}  \\
  \gamma&=arcsin\dfrac{r\cos(\beta-\phi)}{L}=arcsin\dfrac{r\cos(\beta-\phi)}{\sqrt{D^2+r^2+2Dr\sin(\beta-\phi)}}
    \end{align}
$$
!!#009900 !24 由此可知，扇形束卷积反投影算法只需在平行束算法的基础上作适当的加权修正即可得到，实现步骤分为以下3步：!!!
1. **投影函数的修正**   假定在$beta$角下的投影函数是$p_\beta(\gamma)$，若以等角度$\Delta\gamma$采样，则有$\color{red} p_\beta(\gamma)=p_\beta(n\Delta\gamma)$。若探测器通道个数为$N_d$，则$n=-\dfrac{N_d+1}{2}+1 \quad : \quad  \dfrac{N_d+1}{2}-1$。<font color=red> 修正后的投影函数为：</font>$$
  \bbox[#f9de9f,6pt]{
 P_\beta(n\Delta\gamma)=p_\beta(n\Delta\gamma)D\cos(n\Delta\gamma)
 }
 $$
 2. **卷积（滤波）运算** 将<font color=red>修正后的投影函数与滤波函数</font>$h^\prime(\gamma)$作卷积运算，得到：
 $$
\bbox[#f9de9f,6pt]{C_\beta(n\Delta\gamma)=P_\beta(n\Delta\gamma) * h^\prime(n\Delta\gamma),\quad h^\prime(n\Delta\gamma)=\dfrac{(n\Delta\gamma)^2}{2\sin^2(n\Delta\gamma)}h(n\Delta\gamma)}
$$
若采用的是R-L滤波器，则：
$$
h(n\Delta ) = \left\{ {\begin{array}{*{20}{c}}
{\frac{1}{{4{\Delta ^2}}}}&{n = 0}\\
0&{n{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} is{\kern 1pt} {\kern 1pt} {\kern 1pt} even{\kern 1pt} {\kern 1pt} number{\kern 1pt} {\kern 1pt} }\\
{ - \frac{1}{{{\pi ^2}{n^2}{\Delta ^2}}}}&{n{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} is{\kern 1pt} {\kern 1pt} {\kern 1pt} odd{\kern 1pt} {\kern 1pt} number}
\end{array}} \right.{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt}  \Rightarrow {\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} h(n\Delta \gamma ) = \left\{ {\begin{array}{*{20}{c}}
{\frac{1}{{4{{(\Delta \gamma )}^2}}}}&{}\\
0&{{\kern 1pt} {\kern 1pt} }\\
{ - \frac{1}{{{\pi ^2}{n^2}{{(\Delta \gamma )}^2}}}}&{}
\end{array}} \right.
$$
其中$\Delta$为采样周期，在扇形束中采样周期即为$\Delta\gamma$，进一步：
$$
\bbox[#f9de9f,6pt]{
h^\prime(n\Delta\gamma)
\beta()=P_\beta(n\Delta\gamma) * h^\prime(n\Delta\gamma),\quad h^\prime(n\Delta\gamma)=\dfrac{(n\Delta\gamma)^2}{2\sin^2(n\Delta\gamma)}h(n\Delta\gamma)
}
$$
分析$h^\prime(n\Delta\gamma)$在$n=0$时的取值：
$$
\begin{align}
h^\prime(n\Delta\gamma)&=\dfrac{(n\Delta\gamma)^2}{2\sin^2(n\Delta\gamma)}h(n\Delta\gamma)=\dfrac{1}{2}h[\sin(n\Delta\gamma)]  \\
h^\prime(n\Delta\gamma)|_{n=0}&=\dfrac{1}{2}h[\sin(n\Delta\gamma)]|_{n=0}=\dfrac{1}{2}h(0)
\end{align}
$$
于是可得加权后的滤波器卷积核为：
$$
\bbox[#f9de9f,6pt]{
h(n\Delta \gamma ) = \left\{ {\begin{array}{*{20}{c}}
{\frac{1}{{8{{(\Delta \gamma )}^2}}}}&{n = 0}\\
0&{n{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} is{\kern 1pt} {\kern 1pt} {\kern 1pt} even{\kern 1pt} {\kern 1pt} number{\kern 1pt} {\kern 1pt} }\\
{ - \frac{1}{{2{\pi ^2}{{\sin }^2}(n\Delta \gamma )}}}&{n{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} is{\kern 1pt} {\kern 1pt} {\kern 1pt} odd{\kern 1pt} {\kern 1pt} number}
\end{array}} \right.
}
$$

3. **加权反投影** 在第2步中计算出了卷积结果$C_\beta(n\Delta\gamma)$，乘以权重因子$\dfrac{1}{L^2}$，然后反投影：
$$
f(r,\phi) \approx \Delta\beta\sum_{m=1}^M{\dfrac{1}{L^2}C_\beta(n\Delta\gamma)|_{n\Delta\gamma=\gamma=\arcsin\dfrac{r\cos(\beta-\phi)}{\sqrt{D^2+r^2+2Dr\sin(\beta-\phi)}}}}
$$
# 总结
- 平行束射线在滤波的时候，只需要对某个旋转角度$\theta$下，投影得到的数据$t|_{\theta=const}$，进行滤波即可；
- 扇形束射线在滤波的时候，也需要对某个旋转角度$\beta$下，投影得到的数据进行滤波，即对$-\gamma_m \sim \gamma_m|_{\beta=const}$射线的投影数据进行滤波。但由于是**扇形束**，需要对**滤波器卷积核**进行修正$\left( h^\prime(\gamma)=\dfrac{\gamma^2}{2\sin^2\gamma}h(\gamma) \right)$，也需要对$-\gamma_m \sim \gamma_m|_{\beta=const}$**射线投影数据进行修正**$\left( P_\beta (\gamma) = p_\beta (\gamma) D\cos (\gamma) \right)$；
-  经过卷积滤波后的**投影数据**$\dfrac{1}{L^2}C_\beta(\gamma)$就可以用来进行图像重建，当然实际的数据是离散化的，只有有限条射线，覆盖整个截面；
- **图像重建** 我们先给出图像的某个点$M$坐标$(x,y)\rightarrow(r,\phi)$，由此便可以计算出是哪条射线经过了$M$点（用扇束角$\gamma$来标识射线，$\gamma=\arcsin\dfrac{r\cos(\beta-\phi)}{\sqrt{D^2+r^2+2Dr\sin(\beta-\phi)}}|_{\beta=const}$），如果射线不存在，那么利用相邻两条射线的投影数据进行**插值**，将得到的数据**“回填”**到$M$点；
- 更新$\beta$角，重复以上过程。

---

---

# 等距离扇形束（Equal-Angle Fan Beam）CT
下图为等距扇束投影系统的几何结构说明。$D_1D_2$为探测器阵所在位置。$S_oB$为某一射线，与探测器阵相交于点B，探测器阵中点为$Q_s$。为简化推导过程中的数学表达式可设想将$D_1D_2$，平移到正好穿过坐标原点的位置$D^\prime_1D^\prime_2$。它是$D_1D_2$的镜像，称为虚拟探测器，二者关系可由扇形的几何尺寸确定，故**射线的相对位置**也可由$OA$定出（!!#ff0000 **虚拟探测器**!!），线段$OA$的长度是探测器$D^\prime_1D^\prime_2$上的距离$ \color{red}s$。**射线的相对位置**还可以由实际探测器$D_1D_2$的位置$Q_sB$决定。本文利用虚拟探测器来决定射线的位置，因此，射线的投影函数可记为$\color{red}p^f_\beta(s)=p^f(s,\beta)$。

![等距扇束几何结构_1.png](https://docimg6.docs.qq.com/image/AgAABenPl4Ep66Y2HQlO-4DSsyRxUuH5.png?w=1052&h=967)

由几何关系可知，平行束（$t,\theta$）和等距扇形束（$s,\beta$）之间的转换关系为：
$$
\begin{align}
\theta &= \beta+\gamma=\beta+\arctan \dfrac{s}{D}  \\
t        &= s\cos \gamma=s\dfrac{D}{\sqrt{D^2+s^2}}=\dfrac{sD}{\sqrt{D^2+s^2}}
\end{align}
$$
进一步可得：
$$
s=\dfrac{Dt}{\sqrt{D-t^2}}
$$

> 当然我们也可以使用!!#ff0000 真实探测器的参数!!来表示，假设射线源到探测器的距离$S_oQ_s=L$，记$Q_sB=d$，于是，
> $$
t= s\cos \gamma=s\dfrac{L}{\sqrt{L^2+d^2}}
$$
> 还可以进一步整理，由于$s//d$，$\dfrac{s}{d}=\dfrac{D}{L} \Rightarrow s=\dfrac{Dd}{L} $，由几何关系还可以得到$\gamma=\arctan\dfrac{d}{L}$，于是使用实际探测器参数时，平行束（$t,\theta$）和等距扇形束（$d,\beta$）之间的转换关系为：
> $$
\begin{align}
\theta &= \beta+\gamma=\beta+\arctan \dfrac{d}{L}  \\
t        &= s\cos \gamma=s\dfrac{L}{\sqrt{L^2+d^2}}=\dfrac{Dd}{L}\dfrac{L}{\sqrt{L^2+d^2}}=\dfrac{Dd}{\sqrt{L^2+d^2}}
\end{align}
$$

已知平行束重建算法为（**卷积形式**）：
$$
f(x,y) = \int_0^\pi  {\int_{ - {t_m}}^{{t_m}} {p(t,\theta )h(x\cos \theta  + y\sin \theta  - t)dt} d\theta }
$$

其中$h(\cdot)$是滤波器。为了方便推导，换为极坐标形式，令$x=r \cos \phi, y=r\sin \phi$，则$x\cos \theta+y\sin \theta=r\cos(\theta-\phi)$，于是：
$$
f(r,\phi)=\dfrac{1}{2}\int_0^{2\pi}{\int_{ - {t_m}}^{{t_m}} {p(t,\theta )h[r\cos(\theta-\phi)  - t]dt} d\theta}
$$
利用$(t, \theta)\Leftrightarrow(s, \beta)$之间的关系，将平行束的变量替换成扇形束的变量，替换后的微分关系可以利用雅克比行列式表示：
$$
\left| J \right| = \left| {\begin{array}{*{20}{c}}
{\frac{{\partial t}}{{\partial s }}}&{\frac{{\partial t}}{{\partial \beta }}}\\
{\frac{{\partial \theta }}{{\partial s }}}&{\frac{{\partial \theta }}{{\partial \beta }}}
\end{array}} \right| = \left| {\begin{array}{*{20}{c}}
{\dfrac{D^3}{(D^2+s^2)^{3/2}} }&0\\
\frac{{\partial \theta }}{{\partial s }}&1
\end{array}} \right| = \dfrac{D^3}{(D^2+s^2)^{3/2}} 
$$
即$dtd\theta=\dfrac{D^3}{(D^2+s^2)^{3/2}}  ds d\beta$，于是可得：
$$
\begin{align}
f(r,\phi)=& \dfrac{1}{2}\int_{-\arctan(s/D)}^{2\pi -arctan(s/D)}{\int_{ - s_m}^{s_m} {p(\dfrac{sD}{\sqrt{D^2+s^2}},\beta+\arctan \dfrac{s}{D} ) }}\\
&\cdot h[r\cos(\beta+\arctan \dfrac{s}{D}-\phi)  - \dfrac{sD}{\sqrt{D^2+s^2}}]\cdot \dfrac{D^3}{(D^2+s^2)^{3/2}} ds d\beta
\end{align}
$$
考虑到$-\arctan(s/D) \sim 2\pi-\arctan(s/D) $已经覆盖了物体的整个$360^\circ$，所以等效于从0到$2\pi$进行积分，且$p(\dfrac{sD}{\sqrt{D^2+s^2}},\beta+\arctan \dfrac{s}{D} )=p^f(s,\beta)$，因此可得扇形束图像重建算法为：

$$
% \boxed{
% f(r,\phi)=\dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - \gamma_m}^{\gamma_m} {p^f(\gamma,\beta)\cdot
% h[r\cos(\beta+\gamma-\varphi)  - D\sin\gamma]\cdot D\cos\gamma d\gamma} d\beta}
%}
 \bbox[#f9de9f,6pt]{
 f(r,\phi)=\dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - s_m}^{s_m} {p^f(s,\beta)\cdot
 h[r\cos(\beta+\arctan \dfrac{s}{D}-\phi)  - \dfrac{sD}{\sqrt{D^2+s^2}}]\cdot \dfrac{D^3}{(D^2+s^2)^{3/2}} ds }d\beta}
 }  
$$
**但上述公式对$s$的积分不是!!#e60000 卷积!!的形式。**

为了配合卷积形式，根据下图对相关变量重新命名：$A$为重建点，重新定义虚拟探测器上射线位置为$s^\prime$，射线源到重建点的长度为$L$。

![等距离扇形束几何结构_2.png](https://docimg8.docs.qq.com/image/AgAABenPl4Heed3qCl5GS59Wnoywi6P0.png?w=772&h=609)

由相似三角形可得：$\dfrac{s^\prime}{\overline{AP}}=\dfrac{\overline{S_oO}}{\overline{S_oO}+\overline{OP}}=\dfrac{D}{D+\overline{OP}}$，$\overline{OP}=r\sin(\beta-\phi)$，于是
<br/>
$r\cos(\beta-\phi)=\overline{AP}=\dfrac{D+r\sin(\beta-\phi)}{D}s^\prime$，定义：

$$
 \bbox[#f9de9f,6pt]{
 {\color{red}U}={\color{red}\dfrac{D+r\sin(\beta-\phi)}{D}}=\dfrac{D+r\sin\beta\cos\phi-r\cos\beta\sin\phi}{D}=\color{red}\dfrac{D+x\sin\beta-y\cos\beta}{D}
 }
$$

于是，

$$
\begin{align}
  & r\cos(\beta+\arctan \dfrac{s}{D}-\phi)  - \dfrac{sD}{\sqrt{D^2+s^2}}  \\
=& r\cos(\beta-\phi)\cos(\arctan\dfrac{s}{D})-r\sin(\beta-\phi)\sin(\arctan\dfrac{s}{D})- \dfrac{sD}{\sqrt{D^2+s^2}}  \\
= & Us^\prime\dfrac{D}{\sqrt{s^2+D^2}} - r\sin(\beta-\phi)\dfrac{s}{\sqrt{s^2+D^2}} - \dfrac{sD}{\sqrt{D^2+s^2}}   \\
= & \dfrac{s^\prime UD}{\sqrt{s^2+D^2}} - [r\sin(\beta-\phi)+D]\dfrac{s}{\sqrt{s^2+D^2}} \\
\mathop  = \limits^{\color{red}r\sin(\beta-\phi)+D=UD} & \dfrac{s^\prime UD}{\sqrt{s^2+D^2}}-\dfrac{sUD}{\sqrt{s^2+D^2}}
\end{align}
$$

因此，等距离扇形束图像重建算法可以表示为：

$$
 \bbox[#f9de9f,6pt]{
f(r,\phi)=\dfrac{1}{2}\int_{0}^{2\pi}{\int_{ - s_m}^{s_m} {p^f(s,\beta)\cdot
 h[(s^\prime-s)\dfrac{UD}{\sqrt{s^2+D^2}}]\cdot \dfrac{D^3}{(D^2+s^2)^{3/2}} ds }d\beta}
 }
$$

> 斜坡滤波器卷积核的缩放性质：$h(at)=\dfrac{1}{a^2}h(t)$。在等距离扇形束重建公式中，滤波器项为： $$ h\!\left[(s' - s) \cdot \frac{UD}{\sqrt{s^2 + D^2}}\right] $$ 令：
>  - $ t = s' - s $ 
>  - $ a = \dfrac{UD}{\sqrt{s^2 + D^2}} $ 
>  则：
>   $$
 h(at) = \frac{1}{a^2} h(t) = \frac{1}{\left( \dfrac{UD}{\sqrt{s^2 + D^2}} \right)^2} \cdot h(s' - s) = \frac{s^2 + D^2}{U^2 D^2} \cdot h(s' - s) 
 $$
>  **虽然 $a$ 是 $s$ 的函数，但在积分中对每一个 $s$ 值，$a(s)$ 是一个常数（相对于变量 $t=s^\prime−s$ 而言）**，因此缩放性质仍然适用。

代入该性质，

$$
  \begin{align}
f(r,\phi) &= \frac{1}{2}\int_{0}^{2\pi}\int_{-s_m}^{s_m} p^f(s,\beta)\cdot
\left[ \frac{(s^2 + D^2)}{U^2 D^2} h(s' - s) \right]\cdot \frac{D^3}{(D^2+s^2)^{3/2}} \, ds \, d\beta  \\
           &= \frac{1}{2}\int_{0}^{2\pi}\int_{-s_m}^{s_m} p^f(s,\beta)
 h(s' - s)\frac{D}{U^2\sqrt{D^2+s^2}} ds \, d\beta 
\end{align}
$$
由于$U$与$s$无关，进一步整理后可得扇形束的!!#ff0000 **卷积形式**!!滤波反投影算法：

$$
 \bbox[#f9de9f,6pt]{
 \begin{align}
f(r,\phi) &=\int_{0}^{2\pi}
\dfrac{1}{U^2}\left(\int_{-s_m}^{s_m}\left[ \frac{D}{\sqrt{D^2+s^2}}p^f(s,\beta)\right]\cdot\left[
  \frac{1}{2}h(s' - s) \right]ds\right) d\beta  \\
           &=\int_{0}^{2\pi}\dfrac{1}{U^2}\left[ \frac{D}{\sqrt{D^2+s^2}}p^f(s,\beta)\right] * h^\prime(s)d\beta,\quad \color{red} h^\prime(s)=\dfrac{1}{2}h(s)
\end{align}
}
$$
其中$U=\dfrac{D+r\sin(\beta-\phi)}{D}=\dfrac{D+x\sin\beta-y\cos\beta}{D}$，$s$为经过待建点$(r,\phi)$的射线，由几何关系得（图中$s^\prime$即为$s$）：

$$
\begin{align}
\dfrac{s}{\overline{AP}}&=\dfrac{\overline{S_oO}}{\overline{S_oO}+\overline{OP}}=\dfrac{D}{D+\overline{OP}} \\
\overline{OP}&=r\sin(\beta-\phi)  \\
\overline{AP}&=r\cos(\beta-\phi)
\end{align}
$$

于是可得经过待建点$(r,\phi)$的射线$s$为

$$
 \bbox[#f9de9f,6pt]{
 s=\dfrac{D}{D+\overline{OP}}\overline{AP}=\dfrac{D}{D+r\sin(\beta-\phi)}r\cos(\beta-\phi)=\dfrac{r\cos(\beta-\phi)}{U}=\color{red}\dfrac{x\cos\beta+y\sin\beta}{U}
 }
$$

!!#009900 !24 实现步骤分为以下3步：!!!
1. **投影函数的修正**   假定在$beta$角下的投影函数是$p_\beta(s)$，<font color=red> 修正后的投影函数为：</font>
$$
 \bbox[#f9de9f,6pt]{
 P_\beta(s)=p_\beta(s)\frac{D}{\sqrt{D^2+s^2}}
 }
$$

若以等间距$\Delta s$采样，探测器通道个数为$N_d$，则$n=-\dfrac{N_d+1}{2}+1 \quad : \quad  \dfrac{N_d+1}{2}-1$，则<font color=red>采样后的修正投影函数为：</font>
$$
  \bbox[#f9de9f,6pt]{
 P_\beta(n\Delta s)=p_\beta(n\Delta s)\frac{D}{\sqrt{D^2+(n\Delta s)^2}}
 }
 $$
 2. **卷积（滤波）运算** 将<font color=red>修正后的投影函数与滤波函数</font>$h^\prime(s)$作卷积运算，得到：
  $$
\bbox[#f9de9f,6pt]{C_\beta(s)=P_\beta(s) * h^\prime(s),\quad h^\prime(s)=\dfrac{1}{2}h(s)}
$$
同理，离散化后可得：
 $$
\bbox[#f9de9f,6pt]{C_\beta(n\Delta s)=P_\beta(n\Delta s) * h^\prime(n\Delta s)}
$$
若采用的是R-L滤波器，则：
$$
h(n\Delta ) = \left\{ {\begin{array}{*{20}{c}}
{\frac{1}{{4{\Delta ^2}}}}&{n = 0}\\
0&{n{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} is{\kern 1pt} {\kern 1pt} {\kern 1pt} even{\kern 1pt} {\kern 1pt} number{\kern 1pt} {\kern 1pt} }\\
{ - \frac{1}{{{\pi ^2}{n^2}{\Delta ^2}}}}&{n{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} is{\kern 1pt} {\kern 1pt} {\kern 1pt} odd{\kern 1pt} {\kern 1pt} number}
\end{array}} \right.{\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt}  \Rightarrow {\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} h^\prime(n\Delta s ) = \left\{ {\begin{array}{*{20}{c}}
{\frac{1}{{8{{(\Delta s )}^2}}}}&{}\\
0&{{\kern 1pt} {\kern 1pt} }\\
{ - \frac{1}{{{2\pi ^2}{n^2}{{(\Delta s )}^2}}}}&{}
\end{array}} \right.
$$

3. **加权反投影** 在第2步中计算出了卷积结果$C_\beta(n\Delta s)$，乘以权重因子$\dfrac{1}{U^2}$，然后反投影：
$$
f(r,\phi) \approx \Delta\beta\sum_{m=1}^M{\dfrac{1}{U^2}C_\beta(n\Delta s)|_{n\Delta s=s=\dfrac{r\cos(\beta-\phi)}{U}=\dfrac{x\cos\beta+y\sin\beta}{U}}}
$$
这里需要特别说明的是，我们已知的是!!#ff0000 **离散**!!的$C_\beta(n\Delta s)|_{s=-n,\cdots,0,\cdots,n-1,n}$即$C_\beta(-n),\cdots,C_\beta(0),\cdots,C_\beta(n-1),C_\beta(n)$。而通过待重建图像的坐标点$(r,\phi)$或$(x,y)$计算得到的$s=\dfrac{r\cos(\beta-\phi)}{U}=\dfrac{x\cos\beta+y\sin\beta}{U}$却是连续的，因此需要利用离散的数据进行插值运算。































