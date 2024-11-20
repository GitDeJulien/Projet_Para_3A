# Mise en oeuvre d'une méthode de décomposition de domaine de type Schwarz sur maillage cartésien régulier

Ce code à pour but de résoudre l'équation de conduction instationnaire en 2D, $\partial_t u(x,y,t) - D\Delta u = f(x,y,t)$ avec les conditions aux limites $u_{\Gamma_0} = g(x,y,t)$ et $u_{\Gamma_1} = \mathcal{h}(x,y,t)$ en se placant dans le domaine réctangulaire $[0, L_x ] \times [0, L_y ]$ de $\mathbb{R}^2$ avec $\Gamma_0$ le bord regroupant les cotés haut et bas et $\Gamma_1$ le bord regroupant les cotés droite et gauche.

# Compilation and execution

 - Ouvrir un terminal
 - Se placer dans le répertoire 'code_para_3A/'
 -  Compilation :
    - make
 - Execution : 
    - make exec

# Architecture du code

## data.cpp

Lis dans le fichier _data.dat_ toutes les données et les regroupe dans une classe __Data__.
Cette class comport aussi une méthode permettant d'afficher tout les paramètres de simulation

## function.cpp

Initialise une class regroupant toutes les fonctions utiles,

 - Condition initial :  __InitialCondition__
 - Terme source : __SourceFunction__
 - Solution exact si on la connais (pour la validation) : __ExactSolution__
 - Condition de bords : 
    - Droite/Gauche : __BoundaryCondition_h__
    - Haut/Bas : __BoundaryCondition_g__

> [!NOTE]
> Toute ces fonctions sont liées à des __key__, présente dans le fichier _data.dat_, permettant de changer de condition initial, de terme source et de conditions de bord.


# Validation

Afin de valider le code des cas tests ont été implémentés. 

## Cas 1 (Stationnaire)

 - Terme source (key 1) : $f = 2(x-x^2 + y-y^2)$ 
 - Conditions de bord (key 1) : $g=0$ et $\mathcal{h}=0$
 - Solution exacte : $u_{exact} = x(1-x)y(1-y)$

## Cas 2 (Stationnaire)

 - Terme source (key 2) : $f = \sin(x) + \cos(y)$ 
 - Conditions de bord (key 3) : $g=\sin(x) + \cos(y)$ et $\mathcal{h} =\sin(x) + \cos(y)$
 - Solution exacte : $u_{exact} = \sin(x) + \cos(y)$

## Cas 3 (Instationnaire périodique)

 - Terme source (key 3) : $f = e^{(-(x-\frac{L_x}{2})^2)}e^{(-(y-\frac{L_y}{2})^2)}\cos(\frac{\pi}{2}t)$ 
 - Conditions de bord (key 2) : $g=0$ et $\mathcal{h} =1$
 - Pas de solution annalytique trouvée