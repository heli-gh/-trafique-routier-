# Projet TP -- méthode de numérique pour EDO

##                                                                        Pan Ruidong&LI Hengshuo

### Le Sommaire

[TOC]
>Toutes les unités de vitesse ci-desous sont m/s, toutes les unités de distance ci-desous sont m. 
>
>les graphes sont un peu hors de position mais pas très loin, pardonne-nous s'il le vous plaît.

### Introduction
  Dans ce TP, nous allons étudier deux modèles différents de trafic routier et déterminer lequel est le plus réaliste grâce à la simulation. De plus, en modifiant les paramètres du modèle tels que le nombre de voitures, la vitesse maximale autorisée, etc. nous pouvons comprendre comment un bouchon se forme et nous pouvons aussi étudier les conditions dans lesquelles des embouteillages peuvent être efficacement évitée. Les modèles consistent en un système EDO et avec l'aide de cours analyse numérique nous avons étudié la méthode d'Euler, RK et la méthode multipas, nous allons donc utiliser ces trois méthodes pour résoudre les modèles. 

### Un premier modèle d'ordre 1

#### Situation

  n voutures roulant les unes derrière les autres ,$x_i(t)$: la position au cours du temps de la $i^{ème}$ voiture ,  $x_i'(t)$ la vitesse de la $i^{ème}$ voiture . 

$$
x_i' (t) = V (x_i(t), x_ {i+1} (t) − x_i(t)) ,∀i ∈ {0, · · · , n − 1}
$$

  Typiquement, on prendra la fonction V de la forme suivante, où ''d'' est la distance entre chaque voiture, ''$d_{min}$'' est la valeur minimale dans les vleurs de ''d'', "$V_{max}$" est la vitesse maximale autorisée.
$$
V(x,d)=\left\{

            \begin{array}{11}

            0 & si \quad d\le d_{min}\\

            V_{max}\frac{(d-d_{min})^2(5d_{min}-d)^2}{(4d^2_{min})^2}                 &     si \quad d_{min} \ge d \le 3d_{min}\\

            V_{max} & si \quad d\ge d_{min}\\

            \end{array}

            \right.
$$


  Ainsi, nous pouvons avoir les aspects pratiques ci-dessus pour former le système de'EDO sous la forme $x' (t) = f(t, x (t)) \quad x \in \mathbb{R^n} $ en précisant la fonction $f : \mathbb{R^n} \rightarrow \mathbb{R^n}$ 

$ x (t)=(x_0(t),..., x_i(t), x_ {i+1} (t),... ,x_{n-1}(t)) ,∀i ∈ {0, · · · , n − 1}$

$x' (t)=(x'_0(t),..., x'_i(t), x'_ {i+1} (t),... ,x'_{n-1}(t)) ,∀i ∈ {0, · · · , n − 1}$

$x_i' (t) =f_i(t,x(t))= V (x_i(t), x_ {i+1} (t) − x_i(t)) ,∀i ∈ {0, · · · , n − 1} $

#### Théorie

  Nous allons dans la partie théorique pour prouver si notre EDO satisfait avec unique solutions. Nous savons par les conditions données que $ V : \mathbb{R^2} \rightarrow \mathbb{R}$ de la classe $\mathcal{C^1} $, tout les dérivees sont bornées, donc on peut montre que le problème de Cauchy admet une unique solution. 
$$
\begin{eqnarray}
x' = V (x, d)\quad d_{i} = x_{i+1}(t)- x_i(t) \\
\frac{\partial V(x,d)}{\partial x} \stackrel{<}{=} K, \in \mathcal{C^0}\\ 
\frac{\partial V(x,d)}{\partial d} \stackrel{<}{=} K\in \mathcal{C^0}\\
|| V(x_i,d_i)- V(\hat{x}_i,\hat{d}_i)|| &=& |V(x_i,d_i)- V(\hat{x}_i,\hat{d}_i)|\\
&=& \frac{\partial V}{\partial x}(x^*,d^*)(x_i-\hat{x}_i)+\frac{\partial V}{\partial d} (x^*,d^*)(d_i-\hat{d}_i)\\
&\stackrel{<}{=}& K(x_i-\hat{x}_i) +K(d_i-\hat{d}_i)\\
&\stackrel{<}{=}& K||(x_i,d_i)-(\hat{x}_i,\hat{d}_i)||\\

\end{eqnarray}
$$
V est dit lipschitzienne, $x_i' (t) =f_i(t,x(t))= V (x_i(t), x_ {i+1} (t) − x_i(t)) ,∀i ∈ {0, · · · , n − 1}$  donc $f_i(t,x(t))$ est lipschitzienne:    
$$
\begin{eqnarray}
||f_i(t,x) - f_i(t,\hat{x})|| & \stackrel{<}{=} & K||x-\hat{x}||\quad x ,\hat{x}\in \mathbb{R^n}\\
||f(t,x)-f(t,\hat{x})|| & \stackrel{<}{=}& L||x-\hat{x}||
\end{eqnarray}
$$
$ f\in \mathcal{C^1}$ et $f$ est lipschitzienne donc $x' (t) =f(t,x(t)), x(t_0)=x_0$ admet une unique solution.

  Prenons une équation aléatoire et testons-la pour voir si nous pouvons obtenir une solution unique, par exemple : $V_{max}(x) = x + 1 \quad x_0 =0$  
$$
\begin{eqnarray}
\frac{dx}{dt}&=& x+ 1 \\
x &=&Ce^t - 1\quad x_0=0\\
x &=&e^t -1
\end{eqnarray}
$$
donc on a le solution unique $x =e^t -1$

#### Algorithme

Nous avons le schéma à 1 pas $RK$:

$k_{ n,1} = f (t_n , x_n )$                                                  avec   $x_{n+1} = x_n + ∆t(c_1 k_{n,1} + c_3 k_{n,3})$
$k_{n,2} = f (t_n + a\Delta t, x_n + \Delta tak_{n,1} ) $
$k_{n,3} = f (t_n + b∆t, x_n + ∆tbk_{n,2} )$

donc il nous faut vérifier le schéma à 1-pas $RK$ est une méthode de $RK$ en donnant le tableau associé.

  méthode de $RK$ :

$$
\begin{eqnarray}
              x_{n+1} &= &x_n + h_n\phi(t_n ,x_n ,h_n )\\
  				\phi(t,x,h)&=&\sum_{i=1}^{r}c_if(t+\theta_ih,\hat{x}_i)\\
  	              k_{n,i}&=& f(t+\theta_ih,\hat{x}_i)\\
  			 \hat{x}_i&=& x + h\sum_{j=1}^{r}a_{i,j}f(x+\theta_jh,\hat{x}_i)
 \end{eqnarray}
$$
  $r=3 , \theta_1 =0; \theta_2 =a ; \theta_3 =b, h_n=\Delta t $

  donc le schéma à 1-pas $RK$ est une méthode de $RK$ et réprésenter par le tableau suivant:

| $\theta_1=0$ | $0$   | $0$         | $0$   |
| ------------ | ----- | ----------- | ----- |
| $\theta_2=a$ | $a$   | $0$         | $0$   |
| $\theta_3=b$ | $0$   | $b$         | $0$   |
|              | $c_1$ | $c_{2}$=$0$ | $c_3$ |

mais nous devons déteminer le valeur de $a,b,c1,c3$  pour écrire notre algorithme RK en l'ordre maximal .  nous commencons à l'ordre $r=3$. e:matrice initiale; $\theta$ : la matrice diagonale dont la taille est r*r et les éléments diagonaux sont les $\theta_i$;  $c$ le vecteur de $\mathbb{R^r}$ de composantes $c_1,c_2,c_3$

$$
\begin{eqnarray}
 k_{n,i}(h)&=&f(t+\theta_ih,x_n + h\sum_{j=1}^{3}a_{i,j}k_{n,j}(h)) \quad pour\quad i= 1,2,3\\
\phi(t,x,h)&=&\sum_{i=1}^{3}c_ik_{n,i}(h)\\

k_{n,i}(0)& = &f(t,x_n),\\
k_{n,i}'(0)&= &\theta_if'_t(t,x_n)+\sum_{j=1}^3a_{i,j}k_{n,j}(0)f'_x(t,x_n)\\
k_{n,i}''(0)&=& \theta_i^2f''_{t^2}(t,x_n)+2\sum_{j=1}^3\theta _i a_{i,j} k_{n,j}(0)f''_{tx}(t,x) + (\sum^3_{j=1}a_{i,j} k_{n,j}(0))^2f''_{x^2}(t,x_n) + \\&&2\sum^3_{j=1}a_{i,j}k_{n,j}'(0)f'_x(t,x_n)\\

\phi (t,x_n,0)& =& \sum_{i=1}^3c_ik_{n,i}(0)=c^Tef(t,x_n)\\
\phi' (t,x_n,0) &= &\sum_{i=1}^3c_ik'_{n,i}(0)=c^T\theta ef'_t(t,x_n)+ c^TA e f(t,x_n)f'_x(t,x_n)\\
\phi'' (t,x_n,0) &= &\sum_{i=1}^3c_ik''_{n,i}(0)=c^T\theta^2 ef''_{t^2}(t,x_n)+ 2c^T\theta A e f(t,x_n)f''_{tx}(t,x_n)+\\
&&c^T(Ae)^2f''_{x^2}(t,x_n)(f(t,x_n))^2 + 2c^TA\theta e f'_t(t,x_n)f'_x(t,x_n)+2c^TA^2e(f'_x(t,x_n))^2f(t,x_n)\\
f^0(t,x_n)&= &f(t,x_n)\\
f^1(t,x_n)&=&f'_t(t,x_n)+f(t,x_n)f'_x(t,x_n)\\
f^2(t,x_n)&=&f''_{t^2}(t,x_n)+2f(t,x_n)f''_{tx}(t,x_n)+(f(t,x_n))^2f''_{x²}(t,x_n)+\\
   &&f'_t(t,x_n)f'_x(t,x_n)+f(t,x_n)(f'_x(t,x_n))^2
\end{eqnarray}
$$
la condition est obtenue par une fonction (p-1) fois continûment différentialble en utilisant la relation:

$ \frac{\partial^{p-1} }{\partial h ^{p-1} } \phi(t,x_n,h)|_{h=0} = \frac{1}{p}f^{p-1}(t,x_n)$

donc :
$$
Ae= \theta e; c^Te =1; c^TAe=1/2;c^T\theta^2e = 1/3; c^TA\theta e =1/6
$$
$a=1/3,b=2/3,c_1=1/4,c_3=3/4$

Maintenant on teste s'il vérifie la condition de l'ordre 4 :


$$
\left(\begin{array}{ccc}

c_{1} & c_{2} & c_{3}\end{array}\right)\times\left(\begin{array}{ccc}

0 & 0 & 0\\

0 & a^{3} & 0\\

0 & 0 & b^{3}

\end{array}\right)=b^{3}\times c_{3}=\frac{2}{9}\neq\frac{1}{4}\left(\begin{array}{ccc}

c_{1} & c_{2} & c_{3}\end{array}\right)\times\left(\begin{array}{ccc}

0 & 0 & 0\\

0 & a^{3} & 0\\

0 & 0 & b^{3}

\end{array}\right)=b^{3}\times c_{3}=\frac{2}{9}\neq\frac{1}{4}
$$


et il vérifie pas. Donc le schéma à 1-pas $RK$ est de ordre 3 maximale.

#### Simulation

  Pour résoudre ce système d’EDO en machine, on se donnera une condition initiale t.q.
           $x^0_1 < x^0_1 < · · · < x^0_{n−1} $et $ x^0_{i+1} − x^0_i = 20m $pour tout $ i ∈ {0, · · · , n − 2}$.

  Après compléter les codes, nous avons testé pour un temps de simulation T=100s, dt=0.1s et 0.01s avec un nombre de voitures de 10, 40 et 100. Maintenant nous allons analyser les résultats obtenus et faire des commentaires.

##### Simulation du modèle à l'aide de différentes méthodes



<div align="center">
<figcaption> T=100 dt =0.1 pour 10 voitures avec Euler </figcaption>
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_p_0.031S_eula.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_v_eula.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>



La figure à gauche : Positon (m) des voitures en fonction du temps (s)
La figure à droite :  Vitesse (m/s) des voitures en fonction du temps (s) 



<div align="center">
<figcaption> T=100 dt =0.1 pour 10 voitures avec RK </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_p_0.0786S_RK.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_v_RK.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : Positon (m) des voitures en fonction du temps (s)
La figure à droite :  Vitesse (m/s) des voitures en fonction du temps (s) 



<div align="center">
<figcaption> T=100 dt =0.01 pour 40 voitures avec Euler </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.01_40_p_0.8255S_eula.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.01_40_v_eula.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : Positon (m) des voitures en fonction du temps (s)
La figure à droite :  Vitesse (m/s) des voitures en fonction du temps (s)



<div align="center">
<figcaption> T=100 dt =0.01 pour 40 voitures avec RK </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.01_40_p_2.4116S_RK.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.01_40_v_RK.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>




La figure à gauche : Positon (m) des voitures en fonction du temps (s)
La figure à droite :  Vitesse (m/s) des voitures en fonction du temps (s)

  D'après les figures ci-dessus, lorsque nous choisissons la méthode RK pour la simulation du modèle, nous constatons que son modèle est en grande partie similaire à celui obtenu avec la méthode d'Euler, mais nous remarquons également que le temps d'exécution du programme est différent, la méthode RK prenant un peu plus de temps pour s'exécuter. Par exemple pour dt=0.1 avec 10 voitures l’Euler prend 0.031s mais RK prend 0.683s.

##### L'effet de la valeur de dt sur le modèle



<div align="center">
<figcaption> T=100, 40 voitures en utilisant la méthode d’euler </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_v_eula.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.01_40_v_eula.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>




​                                                                                            Vitesse (m/s) des voitures en fonction du temps (s)

La figure à gauche : dt=0.1
La figure à droite :  dt=0.01




<div align="center">
<figcaption> T=100, 100 voitures en utilisant la méthode de RK </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_100_v_RK.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.01_100_v_RK.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>




​                                                                                         Vitesse (m/s) des voitures en fonction du temps (s)

La figure à gauche : dt=0.1
La figure à droite :  dt=0.01

  D'après les figures ci-dessus, lorsque nous choisissons différentes valeurs de dt, pour dt=0,1, dans la figure vitesse-temps correspondante, nous pouvons constater que le changement de vitesse entre les voitures de 90m/s à 50m/s se comporte de manière abrupte, mais pour dt=0,01, la quantité de changement de vitesse entre les voitures est lisse, nous pouvons donc savoir que plus dt est petit, plus la simulation est détaillée et précise, et nous constatons également que lorsque le dt est plus petit, plus le programme prend du temps à s'exécuter. Par exemple pour la méthode Euler avec 40 voitures dt=0.1 prend 0.092s mais dt=0.01 prend 0.826s.

##### L’effet du nombre de voitures sur le modèle

  Ici nous prenons T=100 dt=0.1 avec la méthode d’euler ou de RK (comme les deux méthodes produisent essentiellement le même modèle dans les même conditions)




<div align="center">
<figcaption> avec 10 voitures en utilisant la méthode d'Euler </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_p_0.031S_eula.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_v_eula.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>




La figure à gauche : Positon (m) des voitures en fonction du temps (s)
La figure à droite :  Vitesse (m/s) des voitures en fonction du temps (s)

  Dans la figure gauche, nous pouvons voir que l’image de position-temps ont trois pentes distinctes représentant 90, 50 et 110 m/s. La pente correspondant à 50 m/s est la plus petite, tandis que la pente correspondant à 110 km/h est la plus grande. Nous pouvons également voir que les images de position-temps de chaque voiture sont parallèles, ce qui peut indiquer que les dix voitures ont la prèsque même vitesse. En outre, l'image est continue et linéaire, ce qui signifie que le changement de vitesse est instantané et qu'il n'y a pas de phase de décélération ou d'accélération.

  Dans l'image de la vitesse en fonction du temps, nous pouvons voir une forte diminution de la vitesse de 90m/s à 50m/s et une forte augmentation de la vitesse de 50m/s à 110m/s, ce qui coïncide également avec l'image de la figure gauche.



<div align="center">
<figcaption> avec 10 voitures en utilisant la méthode d'Euler </figcaption>
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_pv_1_eula.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_pv_2_eula.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_pv_3_eula.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_10_pv_4_eula.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>



​                                                                                          Simulation du ralentissement des 10 voitures

  En observant la vitesse des voitures à différents moments, nous pouvons remarquer que la voiture passe de 90 m/s à 47 m/s en une seconde environ, et de 50 m/s à 110 m/s presque instantanément, ce qui est clairement irréaliste.



<div align="center">
<figcaption> avec 40 voitures en utilisant la méthode RK </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_p_0.2586S_RK.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_v_RK.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : Positon (m) des voitures en fonction du temps (s)
La figure à droite :  Vitesse (m/s) des voitures en fonction du temps (s)

  Nous pouvons voir que la figure de position-temps est similaire au précédent, mais nous remarquons également que la différence de distance entre chaque voiture est également plus prononcée, mais autour de 50 secondes, chaque courbe est plus serrée, ce qui indique que la distance entre la première et la dernière voiture est courte.

  La figure de temps-vitesse met en évidence la chute de la vitesse dans le passage de 90m/s à 50m/s, et la forte augmentation de la vitesse dans le passage de 50m/s à 110m/s, or on peut voir que la dernière voiture a pris plus de temps pour passer la zone de 50m/s.



<div align="center">
<figcaption> avec 40 voitures en utilisant la méthode RK </figcaption>
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_pv_1_RK.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_pv_2_RK.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_pv_3_RK.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_40_pv_4_RK.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>



​                                                                                         Simulation du ralentissement des 40 voitures

  Dans ces figures, nous constatons que les temps de décélération et d'accélération sont toujours très courts. Ce modèle simulé suppose donc que la voiture a atteint la nouvelle limite de vitesse par un fort freinage, mais dans la vie réelle, les conducteurs ont besoin d'un certain temps de réaction à la décélération soudaine de la voiture qui les précède. ce que le modèle ne prend évidemment pas en compte.



<div align="center">
<figcaption> avec 100 voitures en utilisant la méthode d'Euler </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_100_p_0.208S_eula.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_100_v_eula.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : Positon (m) des voitures en fonction du temps (s)
La figure à droite :  Vitesse (m/s) des voitures en fonction du temps (s)



<div align="center">
<figcaption> avec 100  voitures en utilisant la méthode d'Euler </figcaption>
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_100_pv_2_eula.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_100_pv_3_eula.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_100_pv_4_eula.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 1-C/0.1_100_pv_5_eula.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>



​                                                                                      Simulation du ralentissement des 100 voitures

  Nous pouvons constater que les observations pour 40 voitures sont encore plus évidentes avec 100 voitures. Et nous pouvons conclure que le ralentissement n’est pas réalisable.



#### Analyse de convergence



<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/eula.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/RKt.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/eula_picture.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/RKp.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>



  Les images de la première rangée sont des captures d'écran des équations analytiques issues de l'analyse de régression de méthode d'Euler et RK3 au terminal, respectivement. 

  Les images de la deuxième rangée sont les images de la relation entre l'erreur et la taille du pas pour Euler et RK, respectivement.

  Par observation, tout d'abord, tant la méthode d'Euler que celle de RK, notre erreur, c'est-à-dire la valeur expérimentale moins la valeur réelle, augmente à mesure que la taille du pas diminue, pour converger progressivement vers zéro. Cela montre que les deux méthodes sont convergentes et que nous pouvons les utiliser en toute confiance.

  Ensuite, nous étudions leur ordre:

  Selon l'équation dérivée de l'ajustement des moindres carrés dans la figure, nous pouvons voir que l'erreur de la méthode d'Euler est proportionnelle à une fois le carré de la taille du pas, tandis que la valeur de l'erreur de notre méthode RK est proportionnelle à trois fois le carré de notre taille de pas, ce qui est cohérent avec la conclusion obtenue dans notre section théorique.

  Nous vérifions cette conclusion en faisant varier la valeur de la taille du pas pour obtenir plusieurs de résultats.



#### Méthode de multipas d'ordre 3

  Nous avons affiné l'algorithme, et voici ses images de vélocité et de position.



<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/multip10.png" width = "50%"  height="70%"  title="méthode de multipas"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/multiv10.png" width = "50%"  height="70%"  title="méthode de multipas"  align=right  />
</div>




​    Les résultats sont largement similaires à la comparaison précédente, et nous pouvons dire que notre code est correct et que la méthode est correcte.

​    Notre prochaine idée est de se concentrer sur l'analyse de son erreur et les propriétés de convergence.




<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/multipast.png" width = "50%"  height="70%"  title="méthode de multipas"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/multipasp.png" width = "50%"  height="70%"  title="méthode de multipas"  align=right  />
</div>




  Nous pouvons voir que notre erreur diminue au fur et à mesure que la taille du pas augmente et converge progressivement vers zéro, comme précédemment, donc cette méthode converge efficacement.

  On voit que le coefficient absolut de son analytique sont 3, c'est à dire qu'il est de ordre 3, parce que nous utilisons la méthode RK qui est l'ordre 3 pour initialiser $y_1$ et $y_2$. Mais quand nous utilisons la méthde d'Euler pour initialiser $y_1$ et $y_2$ le coefficient absolut de son analytique doit être 2,  parce que l'initialisation avec la méthode d'Euler amplifie l'erreur sur $y_1$ et $y_2$ , ce qui conduit à une erreur plus importante dans la méthode multipas, ce qui conduit à un ordre plus faible.



### Un modèle d’ordre 2

  Ce scénario est différent du précédent, et cette fois-ci nous commençons directement par l'accélération pour la limiter. C'est plus réaliste car le conducteur contrôle directement l'accélérateur, ou en d'autres termes, l'accélération.Notre modèle d'accélération correspond à l'équation suivante:

l'accélération:

$x''_i (t) = A (x_i (t), x _{i+1} (t) − x_ i (t), x '_ i (t), x '_{ i+1} (t) − x'_i (t) , ∀i ∈ {0, · · · , n − 1}$

Voici la formule analytique de l'accélération que nous avons obtenue à partir des données:


​				
$$
A(x,d,x',d')=\left\{
       \begin{array}{11}
      12\frac {(V{max}(x)-x')}{V{max}} & si \quad d\ge d{min} & x'\le V{max}(x)& d' \ge 0  \\

      -\frac{(2d{min}-d)x'}{d^2}         &   si \quad d \le 2d{min} & x'\ge 0 & d' \le 0 \\

       -\frac{(d{min}-d)x'}{d^2}         &   si \quad d \le d{min} & x'\ge 0 \\

       12\frac {(V{max}(x)-x')}{V{max}} & si \quad d\ge d{min} & x'\le V{max}(x)& d' \ge 0  \\

      V_{max} & si \quad x'\ge V{max}(x)& d' \ge 0 \\

      \end{array}
      \right.  
$$


où ''d'' est la distance entre chaque voiture, ''$d_{min}$'' est la valeur minimale dans les vleurs de ''d'', "$V_{max}$" est la vitesse maximale autorisée, "$x^{'}$" est la vitesse et "$d{'}$ " est l'importance de la variation de la distance entre chaque voiture. 				


#### Théorique

- Reformuler le system d'EDO d'ordre 2 sous la forme $y'(t)= f(t,y(t))$ 

​          
$$
y(t)= (x(t), x'(t)) \in \mathbb{R^{2n}} x(t)\in \mathbb{R^n},x'(t)\in \mathbb{R^n}\\
f:\mathbb{R^{2n}} \rightarrow \mathbb{R^{2n}}\\
y_i' (t) =f_i(t,y(t))= y_{i+n} ,∀i ∈ {0, · · · , n − 1}\\
y_i'(t) =f_i(t,y(t))=  A (x_i (t), x _{i+1} (t) − x_ i (t), x '_ i (t), x '_{ i+1} (t) − x'_i (t) ,∀i ∈ {n-1, · · · , 2n − 1}\\
$$

#### Simulation
  Notre position initiale est la même que dans le cas précédent, mais nous avons augmenté la vitesse initiale de chaque voiture à 90 m/s.

  Comme dans la partie 1, nous avons testé pour un temps de simulation T=100s, dt=0.1s, 0.01s et 0.002s avec un nombre de voitures de 10, 40 et 100. Maintenant nous allons analyser les résultats obtenus et faire des commentaires.

  En général, le changement de dt n'a pas d'impact significatif sur le modèle (sauf que le temps calcul), mais le modèle varie pour différents nombres de voitures et aussi pour différents méthodes, et nous commenterons des exemples spécifiques par la suite.



<div align="center">
<figcaption> T=100 en utilisant la méthode d'Euler </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-40-p-1.173429S-euler.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.002-100-p-euler.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : Positon (m) des voitures en fonction du temps (s) pour 40 voitures avec dt=0.1
La figure à droite :  Positon (m) des voitures en fonction du temps (s) pour 100 voitures avec dt =0.002

  Par rapport au modèle de premier ordre, l’image de position-temps de ce modèle montre que la distance entre les voitures augmente avec le temps passe. De plus, la pente de cette image change plus lentement, contrairement au modèle de premier ordre qui change très rapidement.




<div align="center">
<figcaption> T=100,  40 voitures en utilisant la méthode RK </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-40-v-RK.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.002-40-v-RK.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : Vitesse (km/h) des voitures en fonction du temps (s) avec dt=0.1
La figure à droite :  Vitesse (km/h) des voitures en fonction du temps (s) avec dt=0.002




<div align="center">
<figcaption> T=100,  40 voitures en utilisant la méthode d'Euler </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-40-v-euler.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.01-40-v-euler.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : Vitesse (km/h) des voitures en fonction du temps (s) avec dt=0.1
La figure à droite :  Vitesse (km/h) des voitures en fonction du temps (s) avec dt=0.01



<img src="https://i.loli.net/2021/05/14/ZTwIfK76rOyAbke.png" style="zoom: 67%;" />




​                                                                                   Vitesse (km/h) des voitures en fonction du temps (s) avec dt=0.002  

  Pour n'importe quel de dt(0.1/0.01/0.002) et quand le nombre de voitures est de 40, la figure de vitesse-temps de la méthode d'Euler et RK montrent que les derniers voitures ont une vitesse nulle dans la région de 50km/h, il est donc évident qu'il y a des embouteillages sous cette situation. Or nous avons également constaté quelques différences entre ces figures obtenues par la méthode d'Euler et celles obtenues par la méthode RK, mais la tendance générale est cohérente. Pour les différentes méthodes, nous fournirons une analyse spécifique après.

   Maintenant nous prenons T=100 dt=0.1 pour 10 voitures avec la méthode d’Euler :




<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-10-pv-1-euler.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-10-pv-2-euler.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-10-pv-4-euler.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-10-pv-5-euler.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>




​                                                                                          Simulation du ralentissement des 10 voitures



Maintenant nous prenons T=100 dt=0.002 pour 100 voitures avec la méthode RK:




<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.002-100-pv-1-RK.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.002-100-pv-2-RK.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.002-100-pv-3-RK.png" width = "50%"  height="70%"  title="RK"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.002-100-pv-4-RK.png" width = "50%"  height="70%"  title="RK"  align=right  />
</div>



​                                                                                        Simulation du ralentissement des 100 voitures

  Nous pouvons clairement observer que par rapport au modèle de premier ordre, dans ce modèle les voitures ont suffisamment de temps pour réduire leur vitesse de 90km/h à 50km/h et accélérer de 50km/h à 110km/h, ce qui rend le modèle plus raisonnable, mais nous pouvons également constater qu'en raison de la lenteur du changement de vitesse dans la zone où la vitesse maximale est de 50km/h, la partie avant des voitures ne ralentissent pas en dessous de 50km/h.

  Si l'on prend 100 voitures, on constate aussi que dans la zone où la vitesse maximale autorisée est de 90 km/h, la vitesse la plus faible des véhicules situés à l'arrière n'est que de 30 km/h environ.

##### Comparaison pour les différentes méthodes

 Dans les observations précédentes, nous avons constaté qu'en utilisant différentes méthodes de simulation, nous pouvons obtenir différents modèles, nous allons maintenant comparer ces deux méthodes.



<div align="center">
<figcaption> T=100,  dt= 0.1 pour 10 voitures   </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-10-v-euler.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.1-10-v-RK.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>




La figure à gauche : en utilisant la méthode d'Euler
La figure à droite :  en utilisant la méthode RK




<div align="center">
<figcaption> T=100,  dt= 0.01 pour 100 voitures   </figcaption>
<img 
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.01-100-v-euler.png" width = "50%"  height="70%"  title="RK"  align="left"/>
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/ddd/test code 2-a/0.01-100-v-RK.png" width = "50%"  height="60%"  title="RK"  align="right"  />
</div>



La figure à gauche : en utilisant la méthode d'Euler
La figure à droite :  en utilisant la méthode RK

  Par rapport à la méthode RK, pour la méthode d’Euler, il y a plus de voitures roulant à des vitesses plus faibles entre 50s et 70s, ce qui rend la vitesse globale des voitures plus faible et le passage plus lent (en particulier pour 10 voitures). En revanche, pour la méthode RK la vitesse minimale est relativement élevée, moins de voitures circulent à la vitesse minimale ce qui rend la vitesse globale plus élevée, le passage plus rapide et moins susceptible de provoquer des embouteillages. Nous pouvons donc conclure que la méthode RK est meilleure. Et dans l'étude suivante pour l’installation du ralentisseur nous utiliserons la  méthode RK.

#### Ralentisseur

  Nous ajoutons un ralentisseur avant la zone de 50m/s, nous allons comparer l'impact du différent type de ralentisseur.

  Nous utilisons l'image de densité d'analyse qui a été fourni et qui comporte la coordonnée x pour le temps, la coordonnée y pour le voiture or nous ne considérons pas la dernière voiture, et la coordonnée z est l'inverse de la distance entre la voiture arrière et la voiture avant, donc plus grande la valeur de z, plus petite la distance à la voiture avant.

  Nous choisissons trois types de ralentisseurs: le premier à vitesse constante 70,  le deuxième la vitesse maximale doit être conforme à cette équation $y=90-0.02x$, et le troisième à vitesse constante 60.

  Nous allons étudier deux points principaux : 1) la variation de la distance entre les voitures due à la décélération ; 2) la distance de la dernière voiture qui passe pendant ce temps.

##### Variation de la distance entre les voitures

  La raison pour laquelle nous nous penchons sur ce problème est que les changements soudains de vitesse entraînent des changements soudains de distance et donc des accidents  de la circulation. 

 Ici nous choisissons seulement dix voitures et observons le changement de vitesse et la distance entre les voitures.




<div align="center">
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/90v10.png" width = "25%"  height="70%"  title="90"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/70v10.png" width = "25%"  height="70%"  title="70"  align=right />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/60v10.png" width = "25%"  height="70%"  title="60"  align=right />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/1sv10.png" width = "25%"  height="70%"  title="y=90-0.02x"  align=right  />
</div>




  Notre premier graphique à gauche est sans ralentisseur, on peut voir que le changement de vitesse des voitures n'est pas très synchronisé, tandis que pour le deuxième graphique qui est avec le deuxième type de ralentisseur la vitesse pour chaque voiture est très synchronisée. Quant aux troisième et quatrième graphiques, respectivement avec la vitesse constante de 60 et 70 comme ralentisseur, nous constatons que   choisissant la vitesse constante de 70 comme ralentisseur est meilleur, la vitesse de chaque voiture est fondamentalement synchrone, mais celle de 60 est légèrement désynchronisée.




<div align="center">
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/90d10.png" width = "50%"  height="70%"  title="90"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/70d10.png" width = "50%"  height="70%"  title="70"  align=right />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/60d10.png" width = "50%"  height="70%"  title="60"  align=right />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/1sd10.png" width = "50%"  height="70%"  title="y=90-0.02x"  align=right  />
</div>



  Notre premier graphique est sans ralentisseur, et nous pouvons voir que la distance entre les voitures change très radicalement pendant la décélération soudaine, tandis que le deuxième graphique est avec une vitesse constante de 70 comme ralentisseur, et la distance entre les voitures ne diminue pas brusquement comme ralentisseur, et de même pour le quatrième graphique est avec une vitesse constante de 60 comme ralentisseur, et la distance entre les voitures change relativement lentement. Mais le troisième graphique est le plus performant, dans lequel nous utilisons le deuxième type de ralentisseur, la distance entre les voitures ne diminue pas brusquement et la distance finale entre les voitures est plus ou moins la même, restant à 44,5m.

 Si l'on souhaite que le changement de vitesse soit aussi régulier que possible et que la distance entre les voitures soit aussi exempte de diminutions brutales que possible, il est clair qu'il est préférable d'utiliser comme ralentisseur la vitesse maximale satisfaisant à l'équation : $y=90-0.02x$ (c'est-à-dire le deuxième type de ralentisseur), mais dans la vie réelle, il est évidemment difficile de limiter la vitesse maximale qui est sous une équation linéaire, et nous choisissons une constante vitesse de 70 ou 60 comme ralentisseur, car la variation de la vitesse de voiture a évolué de manière régulière et la distance minimale de l'atelier est estimée à 18 m.

##### Distance du dernier véhicule à passer
  Nous allons explorer l'effet du ralentissement sur la dernière voiture, et nous observons la distance parcourue par la dernière voiture .
  Vous trouverez ci-dessous les figures de la distance parcourue par la dernière voiture en fonction du nombre total de voitures, la coordonnée horizontale étant le nombre total de voitures et la coordonnée verticale étant la distance parcourue par la dernière voiture.



<div align="center">
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/90lastcar.png" width = "50%"  height="70%"  title="90"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/70lastcar.png" width = "50%"  height="70%"  title="70"  align=right />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/60lastcar.png" width = "50%"  height="70%"  title="60"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/1slastcar.png" width = "50%"  height="70%"  title="y=90-0.02x"  align=right  />
</div>



  Nous savons déjà que pour 40 voitures sans ralentisseurs il provoque des embouteillages et parce que la vitesse des voitures à l'arrière est réduite à zéro, mais si nous utilisons des ralentisseurs, cette situation n'existera pas du tout.

  Nous présentons d'abord ces figures, la figure 1 est le cas sans ralentisseur, la figure 2 est le cas avec le ralentisseur à une vitesse constante de 70, la figure 3 est le cas avec le ralentisseur à une vitesse constante de 60, et la figure 4 est le cas avec le ralentisseur à une vitesse maximale satisfaisant $y=90-0.02x$.

  En comparant ces figures, nous pouvons voir que nous avons résolu avec succès le problème de embouteillage en utilisant des ralentisseurs.  

  Dans le cas avec le ralentisseur à une vitesse constante de 70 et à une vitesse maximale satisfaisant $y=90-0.02x$ la dernière voiture parcourt la même distance, probablement parce que la vitesse moyenne des voitures dans la zone de ralentisseur est la même. En comparant la distance parcourue par la dernière voiture avec un ralentisseur à une vitesse constante de 70 et de 60, on constate que plus la vitesse est élevée dans la zone du ralentisseur, plus la distance parcourue par la dernière voiture est long. Nous devons donc choisir un ralentisseur avec une vitesse moyenne plus élevée dans la zone du ralentisseur.  

##### Cas optimale

  Ensuite, nous allons choisir la vitesse optimale pour le ralentisseur. Parce que nous ne pouvons pas installer un ralentisseur à une vitesse maximale satisfaisant $y=90-0.02x$, donc nous cherchons des solutions optimales plus grandes que la vitesse constante 70.

  Nous choisissons des vitesses de 85, 80, 75 et 70 m/s pour la zone du ralentisseur, puis nous comparons les fluctuations de la distance entre chaque voiture, observer s'il y a un embouteillage et la distance parcourue par la dernière voiture.

 Tout d'abord, nous testons s'il y a des embouteillages en utilisant différents types de ralentisseur. Selon l'image, il n'y a pas d'embouteillage. Nous ne montrons pas les images ici.

  Ensuite, nous comparons la variation de la distance entre chaque voiture, les images suivantes sont les ralentisseurs à une vitesse constante de 70, 75, 80 et 85.

  

<div align="center">
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/70d20.png" width = "50%"  height="70%"  title="90"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/75d20.png" width = "50%"  height="70%"  title="70"  align=right />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/80d20.png" width = "50%"  height="70%"  title="60"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/85d20.png" width = "50%"  height="70%"  title="y=90-0.02x"  align=right  />
</div>



  Nous pouvons voir que pour le ralentisseur à une vitesse constante de 80 et 85, la distance entre les voitures varie drastiquement et la distance minimale peut atteindre 15m, et cela peut affecter 15 voitures, donc c'est très dangereux.

  Mais pour le ralentisseur à une vitesse constante de 75, la variation de la distance entre les voitures est relativement faible et la distance minimale est d'environ 18 m, ce qui n'affecte que 7 véhicules, et est donc acceptable.

  Ensuite nous comparons la distance parcourue par la dernière voiture, les figures sont dans le même ordre que ci-dessous :

  

<div align="center">
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/70lastcar.png" width = "50%"  height="70%"  title="90"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/75lastcar.png" width = "50%"  height="70%"  title="70"  align=right />
</div>
<div align="center">
<img
src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/80lastcar.png" width = "50%"  height="70%"  title="60"  align=left />
<img src="/home/heli/Documents/GM6/analysenumérique/EDOTP/Codes/Codes/85lastcar.png" width = "50%"  height="70%"  title="y=90-0.02x"  align=right  />
</div>




   Nous avons précédemment conclu que plus la vitesse moyenne des voitures dans le zoom du ralentisseur est élevée, plus la distance parcourue par la dernière voiture est longe, et évidemment les figures ci-dessus vérifie à nouveau cette conclusion. 

  Or nous constatons que les trois derniers figures ne diffèrent pas beaucoup. Donc tout bien considéré, nous pensons que le meilleur choix de ralentisseur devrait être à une vitesse constante proche de 75. Dans ce cas, le problème du embouteillage est résolu, il y a moins de risques d'accidents de la route et la distance parcourue par la dernière voiture est assez longe.


### Conclusion

  Le premier modèle que nous étudions construit par contrôler la vitesse.  Une fois que nous avons un modèle concret,  nous vérifions que la solution unique est satisfaite sous les conditions de Cauchy en utilisant les théorèmes que nous avons appris en classe.  Après cela, nous avons completé les algorithmes à simuler sur l'ordinateur, et nous avons choisi trois méthodes pour étudier: la méthode d'Euler, le RK, et le multipas. Ensuite nous avons effectué une analyse de régression pour chacune des trois méthodes, ainsi qu'une comparaison de la relation entre l'erreur propre et la taille du pas(dt). Les trois méthodes ont été utilisées pour classer et discuter les cas avec différents nombres de voitures et différents taille du pas. Finalement nous concluons que le modèle n'est pas réaliste.

   Ensuite, nous contrôlons l'accélération pour construire le deuxième modèle, nous construisons toujours le modèle d'abord en conjonction avec la réalité, suivi d'une analyse pour différentes voitures et la taille du pas. Et nous rencontrons des embouteillages, nous explorons ensuite le rôle que jouent les ralentisseurs dans les embouteillages et nous trouvons la limite de vitesse optimale pour le ralentisseur en fonction de la distance entre chaque voiture et de la distance parcourue  par la dernière voiture.

​    Grâce à ce TP, nous avons appris plus intuitivement à résoudre des équations différentielles en utilisant différentes méthodes et nous avons appris à utiliser les ordinateurs pour résoudre des équations différentielles afin d'analyser des problèmes.