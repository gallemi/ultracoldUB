
--------------------------
Resumen de comandos de git
--------------------------


- La primera vez "clonar" el repositorio, usando por ejemplo el 
clone or download de la web de github

- Al unzippear creara un directorio ~/git/ultracoldub este directorio 
contiene el repositorio y toda la informacion sobre ramas, ficheros


- Entrar en ~/git/ultracoldub

- git status

On branch master
Your branch is up-to-date with 'origin/master'.

Untracked files:
  (use "git add <file>..." to include in what will be committed)

	.resumencomandosgit.md.swp

nothing added to commit but untracked files present (use "git add" to track)

> esto significa que he creado un nuevo fichero que git no conoce

- git add resumentcomandosgit.md

ahora git ya lo ha añadido a mi repositorio local, dentro de la rama principal "origin/master"

- git commit  -m "resumen comandos git"

pone al dia el contenido de los ficheros que he modificado dentro de mi repositorio

- git push 

sube al repositorio remoto 

