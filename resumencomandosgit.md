
--------------------------
Resumen de comandos de git
--------------------------


- La primera vez "clonar" el repositorio, usando por ejemplo el 
clone or download de la web de github

- Al unzippear creara un directorio ~/git/ultracoldub este directorio 
contiene el repositorio y toda la informacion sobre ramas, ficheros


- Entrar en ~/git/ultracoldub

- <b>git status</b>

On branch master
Your branch is up-to-date with 'origin/master'.

Untracked files:
  (use "git add <file>..." to include in what will be committed)

	.resumencomandosgit.md.swp

nothing added to commit but untracked files present (use "git add" to track)

> esto significa que he creado un nuevo fichero que git no conoce

- <b>git add</b> resumentcomandosgit.md

ahora git ya lo ha añadido a mi repositorio local, dentro de la rama principal "origin/master"

- <b>git commit</b>  -m "resumen comandos git"

pone al dia el contenido de los ficheros que he modificado dentro de mi repositorio

- <b>git push </b>

sube al repositorio remoto 

---------------------------------
tutoriales
---------------------------------
carla marin dixit, 

Mira, este (acabo de descubrirlo) está bastante bien para familiarizaros con los comandos en la terminal, es básico pero incluye también branches y merges:
https://try.github.io/levels/1/challenges/1

Y para las opciones que tienes desde la web de github su propio tutorial que comentamos ayer también está muy bien:
https://guides.github.com/activities/hello-world/
