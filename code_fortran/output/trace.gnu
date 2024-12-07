# permet de définir les couleurs d'affichage
set palette defined ( 0 '#F7FBFF',\
                      1 '#DEEBF7',\
                      2 '#C6DBEF',\
                      3 '#9ECAE1',\
                      4 '#6BAED6',\
                      5 '#4292C6',\
                      6 '#2171B5',\
                      7 '#084594' )

# pour faire des png. Commenter pour un affichage classique.
#set terminal png

# forcer l'étendue de la colorbar. A changer selon les cas.
#set cbrange [-0.2:1.2]

# force la mise à l'échelle des axes.
set size ratio -1

# si affichage 3D, force l'étendue de l'axe des z. A changer selon les cas.
#set zrange [-0.2:1.2]
# si affichage 3D, fixe le point de vue. A changer selon les cas.
set view 98,140

#pour faire un gif
set term gif animate 
set output "sol.gif"

do for [i = 1:100] {
    t=(i-1)*1.0
    set title "t = ".sprintf("%f", t)." s"." (i = ".sprintf("%d", i).")"
    show title
    # fait un affichage 2D en couleur
    #plot "./sol.".i.".dat" u 1:2:3 palette with image
    ## fait un affichage en 3D et en couleur
    ## Pour le gif
    splot "./sol.".i.".dat" u 1:2:3 palette title "solution approchee", "./sol.".i.".dat" u 1:2:4 title "solution exact"
}

unset term
unset view