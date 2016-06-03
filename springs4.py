#
# Springs v.4
#
#  Version operativa.
#  Funciona con python 2.7
#
#  Ricardo 3-Junio-2016/ 11:00
#          
#
# Dibuja:
#   un resorte vertical
#   un suelo de madera
#   un protipo de regla
#   grafica con variaciones de energia en tiempo real
#   grafica con la aceleracion, la velocidad y la posicion
#
# FALTA:
#
#  Todo el toston de los widgets que esta en fase de desarrollo en la version 5.
#  Hacer la regla no tan cutre
#  Emplear el uso de frames
#  Hay demaseadas cosas 'hardwired' se tendrian que poner en terminos de variables
#  Revisar que las cosas esten en la posicion que se dice.
#
# PREVISION PARA LA VERSION 5
#
#  Que se puedan entrar los datos via widgets
#  La regla con tus etiquetas estara en un frame
#
# POSIBLE CHORRADA A PONER:
#
#  Que en las graficas la escala horizontal vaya corriendo y no sea fijo el cero
#  si no que la grafica se vea como una 'ventana' de tiempo de anchura fija
#  que va corriendo.
#
##########################################################
### Carga las librerias que se emplearan en el programa
##########################################################

from __future__   import division # Prevee problemas con la divison entera del 2.7
from visual       import *        # Carga la libreria VPython
from visual.graph import *        # Carga la parte de graficas del vpython

#...................
# Datos Iniciales
#...................

m    = 0.5                  # Valor de  la masa en         (kg)
R0   = vector(0., 4., 0.)   # Posicion  inicial de la masa (m)
Vel0 = vector(0., 1., 0.)   # Velocidad inicial masa       (m/s)
L0   = 0.                   # Longitud  natural del muelle (m)
k    = 3.                   # Constante recuperadora del muelle
g    = -9.8                 # Valor de la constante de gravedad.


deltat = 0.001              # Paso de tiempo
Tini   = 0.                 # Tiempo inicial
Tmax   = 20.                # Maximo valor del tiempo

#....................................................................
# (A partir de aqui  es logica de programa... no se debiera toquitear
#  salvo para mejora del algoritmo)
#....................................................................


Ventana = display(title="Oscilador harmonico vertical",
                  x=20, y=0, width=420,height=800,)

Ventana.backgroud = color.black     # Fondo del dibujo en color negro
Ventana.center    = (0., 2., 0.)    # Vision centrada en el punto.


Energias   = gdisplay(x=440,y=0,width=800,height=400,
                    xmin=Tini,xmax=Tmax,xtitle="t (s)",
                    ymin=0.,ymax=20.,ytitle="K (J), Ugrav",
                    title="Energies vs temps")

Posiciones = gdisplay(x=440,y=430,width=800,height=400,
                    xmin=Tini,xmax=Tmax, xtitle="t (s)", 
                    ymin=-10.,ymax=10.,ytitle="y (m), vel (m/s), a (m/s**2)",
                    title="Posiciones")


ECinetica  = gcurve(gdisplay=Energias)
EGravitat  = gcurve(gdisplay=Energias)
EElastica  = gcurve(gdisplay=Energias)

YPosition  = gcurve(gdisplay=Posiciones)
Velocity   = gcurve(gdisplay=Posiciones)
Aceleracion= gcurve(gdisplay=Posiciones)


label(display=Energias.display, pos=(0.6*Tmax,0.80*20),  color=color.red,   text="Elastica",box=0,xoffset=1 )
label(display=Energias.display, pos=(0.6*Tmax,0.65*20),  color=color.cyan,  text="Gravitatoria",box=0,xoffset=1)
label(display=Energias.display, pos=(0.6*Tmax,0.50*20),  color=color.yellow,text="Cinetica",box=0,xoffset=1 )

label(display=Posiciones.display, pos=(0.6*Tmax,0.80*10),color=color.red,   text="Posicion",box=0,xoffset=1 )
label(display=Posiciones.display, pos=(0.6*Tmax,0.65*10),color=color.cyan,  text="Velocidad",box=0,xoffset=1)
label(display=Posiciones.display, pos=(0.6*Tmax,0.50*10),color=color.yellow,text="Aceleracion",box=0,xoffset=1 )



# Definicion de los objetos


# MUELLE:

LongInicialMuelle = R0.y-L0                   # Longitud inicial del muelle.

Muelle = helix(axis = (0, 1, 0),
               length = LongInicialMuelle,
               radius = 0.1,
               thickness = 0.05,
               color = color.green,
               ncoils=5)

# FLECHA:

Flecha = arrow(shaftwidth=0.001,
               headlenght=4.,
               headwidth=0.006,
               length=0.5,
               axis=(-0.0001,0.,0.),
               color=color.cyan)

#SUELO y PARED

Suelo = box(pos=vector(0,-0.1,0),  size=(1.,0.1,1.),material=materials.wood)


# RULER CUTRE PERO DE MOMENTO VA. No es mas que un 'prototipo'

def ReglaY(miny,mayu,pasoy):
    Ruler = box(pos=vector(-0.75,2.5,-0.05),size=(0.5,5.3,0.1),material=materials.plastic)
    y=miny
    paso2 = pasoy/5.
    
    tics    = [pasoy*i  for i in arange(0.,mayu+pasoy,pasoy)]
    subtics = [paso2*i  for i in range(1,5)]

    for y in tics:
        label(display=Ventana,
              pos=vector(-1,y,-0.5),
              color=color.red,
              text=str(y)+"--",
              box=0,
              opacity=0,
              xoffset=1 )
        if y == mayu:
            break
        for ys in subtics:
            label(display=Ventana,
            pos=vector(-0.58,y+ys,-0.5),
            color=color.blue,
            text="-",
            box=0,
            opacity=0)
                
 
        
            


# MASA:

masa = sphere(radius = 0.2,
              color = color.yellow)

masa.m   = m                    # Valor de  la Massa
masa.pos = R0-vector(0.,L0,0.)  # Posicion  inicial de la masa
masa.v   = Vel0                 # Velocidad inicial de la masa

Flecha.pos = masa.pos

ReglaY(0.,5.,1.)

# TIEMPOS:

t        = Tini                 # Tiempo inicial.

#...........................................
# BUCLE DE CALCULO SOBRE LA VARIABLE tiempo
#...........................................

while t < Tmax :                                     # Bucle sobre tiempos hasta alcanzar un Tmax
 
    rate(1000)                                       # Frenazo al PC para que los dibujos se vean de forma suave.
    
    EstiramientoMuelle = (Muelle.length - L0 -R0.y)  # Valor de lo que se estira el muelle.

    Fmuelle  = -k * EstiramientoMuelle      # Fuerza del muelle
    Fgrav    = masa.m*g                     # Fuerza gravitatoria
    Fres     = vector(0.,Fmuelle+Fgrav, 0.) # Fuerza resultante
    
    masa.v   = masa.v   + (Fres/masa.m * deltat) # Nueva velocidad (2a. Ley Newton)
    masa.pos = masa.pos + masa.v * deltat        # Nueva posicion.
    Muelle.length = masa.pos.y
    Flecha.pos    = masa.pos
    
    Ecin  = 0.5 * masa.m * mag(masa.v)**2          # Energia Cinetica
    Eelas = 0.5 *   k    * EstiramientoMuelle**2   # Energia Potencial elastica
    Egrav = m * (-g) * masa.pos.y                  # Energia Potencial gravitatoria respecto al suelo
    
    ECinetica.plot(pos=(t,Ecin),       color=color.yellow)  # Ploteo de la energia cinetica
    EElastica.plot(pos=(t,Eelas),      color=color.red)     # Ploteo de la energia elastica
    EGravitat.plot(pos=(t,Egrav),      color=color.cyan)    # Ploteo de la energia gravitatoria
    
    YPosition.plot(  pos=(t,masa.pos.y),color=color.red)    # Ploteo de la posicion
    Velocity.plot(   pos=(t,masa.v.y),  color=color.cyan)   # Ploteo de la velocidad
    Aceleracion.plot(pos=(t,Fres.y/m),  color=color.yellow) # Ploteo de la acelaracion
       
    t = t + deltat #Actualizacion de tiempos
    
        
#.........................................

print "Calculo finalizado"
    
