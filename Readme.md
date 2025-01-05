# Mise en oeuvre d'une méthode de décomposition de domaine de type Schwarz sur maillage cartésien régulier

Ce code à pour but de résoudre l'équation de conduction instationnaire en 2D, $\partial_t u(x,y,t) - D\Delta u = f(x,y,t)$ avec les conditions aux limites $u_{\Gamma_0} = g(x,y,t)$ et $u_{\Gamma_1} = \mathcal{h}(x,y,t)$ en se placant dans le domaine réctangulaire $[0, L_x ] \times [0, L_y ]$ de $\mathbb{R}^2$ avec $\Gamma_0$ le bord regroupant les cotés haut et bas et $\Gamma_1$ le bord regroupant les cotés droite et gauche.

# Compilation et exécution

 - Ouvrir un terminal
 - Se placer dans le répertoire 'code_fortran/' là où le __Makefile__ est présent
 -  Compilation :
    - __make__ ou __make debug__ (pour une compilation en mode _debug version_)
    - __make release__ (pour une compilation en mode _release version_)
 - Execution : 
    - __make exe_debug__ (exécution _debug_)
    - __make exe_release__ (exécution _release_)

# Architecture du code

## sub/mod_toml_parser.f90

Ce module permet de lire et d'afficher les données présentes dans un fichier _.toml_

## sub/mod_precision.f90

Ce module définie le paramètre de précision $pr=8$ (double précision) ainsi que le nombre $\pi$

## sub/mod_charge.f90

Module incluant une subroutine permettant une distribution équitable de la charge et une autre subroutine permettant d'ajouter un recouvrement (overlap) dans le cas d'un maillage structuré.

## src/mod_data.f90

Ce module permet de lire les données dans le fichier _data.toml_ (path = ./data/data.toml) et les regrouper dans un type dérivé

### Contenu

- _type_ __DataType__ : regroupe les données
- _subroutine_ __config_data__ : lis les données dans _data.toml_

## src/mod_functions.f90

Ce module regroupe les fonctions utiles à la définition du problème

### Contenu

 - _function_ __InitialCondition__ : Définie la condition initiale
 - _function_ __ExactSolution__ : Définie la solution exacte si celle-ci ce calcul, sinon renvoie $0.0$
 - _function_ __SourceTerme__ : Définie la fonction source
 - _function_ __BC_Left__ : Définie la fonction de bord gauche
 - _function_ __BC_Right__ : Définie la fonction de bord droit
 - _function_ __BC_Up__ : Définie la fonction de bord haut
 - _function_ __BC_Down__ : Définie la fonction de bord bas

## src/mod_scheme.f90

Ce module permet détablir le produit matrice vecteur qui nous intéresse ainsi que le terme source. Comme on est en 2D nous avons choisit de définir des vecteurs avec une taille correspondant au nombre de points de calcul présent dans chaque sous-domaine. Ce module permet aussi l'initialisation de la solution, le calcul de la solution exacte à chaque pas de temps et l'envoi des message pour la communication MPI. La stratégie de comunication est décrite dans le rapport.

### Contenu

 - _function_ __Lap_MatVectProduct__ : Définie le produit matrice vecteur du laplacien avec résolution Euler Implicite.
 - _function_ __SourceTerme__ : Définie le terme source
 - _subroutine_ __InitSol__ : Définie la solution initiale approchée et exacte si elle est définie.
 - _subroutine_ __ExactSolFunct__ : Définie la solution exacte à chaque temps
 - _subroutine_ __SendMessage__ : Permet l'envoi de messages en fonction des deux conditions au limites définies entre processeur (Dirichlet et Robin)

## src/linear_algebra.f90

Ce module définie un solveur BiCGStab (Bi Conjugate Gradient Stabilized) pour résoudre le problème implicite. Nous avons choisit ce solveur car avec les conditions aux limites de Robin entre les processeurs, la matrice du laplacien n'est pas forcément symétrique.

### Contenu

 - _subroutine_ __Lap_BiCGStab__ : Résout le système AX = b en utilisant le produit matrice vecteur du laplacien définie dans 'mod_scheme.f90'

## src/time_advance.f90

Ce module permet de faire une itération en temps 

### Contenu

 - _subroutine_ __Advance__ : Construit d'abord le terme source puis appelle le BiCGStab pour la résolution.


## src/mod_save_output.f90

Ce module permet de la sauvegarde dans des fichier _.dat_ ou _.vtk_ des données calculées

### Contenu

 - _subroutine_ __SaveSol__ : Sauvegarde la solution approchée
 - _subroutine_ __SaveSolExact__ : Sauvegarde la solution exacte si elle est calculée
 - _subroutine_ __SaveErr__ : Sauvegarde l'erreur en norme 2 entre la solution exacte et approchée 
 - _subroutine_ __SaveTime__ : Sauvegarde le temps d'exécution

## src/main.f90

Le programme avec la boucle en temps


# Validation

Afin de valider le code des cas tests ont été implémentés. Dans le fichier _data.toml_ il suffit de changer le code de __cas__ pour faire le calcul sur les cas tests présentés ci-dessous :

## Cas 1 (Stationnaire)

 - Terme source : $f = 2(x-x^2 + y-y^2)$ 
 - Conditions de bord : $g=0$ et $\mathcal{h}=0$
 - Solution exacte : $u_{exact} = x(1-x)y(1-y)$
 - Condition initiale : $u_{0} =$ 20*SolExact(X,0.0)

## Cas 2 (Stationnaire)

 - Terme source : $f = \sin(x) + \cos(y)$ 
 - Conditions de bord : $g=\sin(x) + \cos(y)$ et $\mathcal{h} =\sin(x) + \cos(y)$
 - Solution exacte : $u_{exact} = \sin(x) + \cos(y)$
 - Condition initiale : $u_{0} = \cos(x)$

## Cas 3 (Instationnaire)

  - Terme source : $f = \sin(x) + \cos(y)$ 
 - Conditions de bord : Dirichlet homogène
 - Solution exacte : $u_{exact} = \sin(\pi x)\sin(\pi y)\exp(-t)$
 - Condition initiale : $u_{0} = \sin(\pi x)\sin(\pi y)$

## Cas 4 (Instationnaire périodique)

 - Terme source : $f = e^{(-(x-\frac{L_x}{2})^2)}e^{(-(y-\frac{L_y}{2})^2)}\cos(\frac{\pi}{2}t)$ 
 - Conditions de bord : $g=0$ et $\mathcal{h} =1$
 - Pas de solution annalytique trouvée
 - Condition initiale : $u_{0} = \exp(-(x-Lx/2)^2)\exp(-(y-Ly/2)^2)$


# Visualisation

Afin de visualiser la solution calculé et la solution exacte il faut :
 - se rendre dans le répertoir _output/_ : __cd output__
 - exécuter le script 'trace_para.sh' : __./trace_para.sh [nombre itération temps]__ (si celui-ci n'a pas été préalablement mis en exécutable taper la commande __chmod +x trace_para.sh__)
 - exécuter le scripte 'trace.gnu' : __gnuplot trace.gnu__ (dans le fichier 'trace.gnu' dans la boucle sur i faire correspondre la borne de fin avec le nombre d'itération en temps)
 - un fichier _sol.gif_ a été créé 

