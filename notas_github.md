### Pasos para actualizar el repositorio local y hacer un pull request con los cambios:

- [Sincronizar el repositorio local (fork) con el original](https://help.github.com/articles/syncing-a-fork/):

      $ git fetch upstream 

      $ git checkout master

      $ git merge upstream/master

- [Actualizar mi repositorio con el local](https://help.github.com/articles/pushing-to-a-remote/):

      $ git push origin master

- Para actualizar el repositorio con los cambios locales:

      $ git status

      $ git add

      $ git commit -m "mensaje"

      $ git push origin master

      $ git status  # para verificar

- El pull request se puede hacer desde la web.


### Sobre branches:

- [Crear una nueva rama localmente (el fork ha de estar actualizado)](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches):

      $ git checkout -b [name_of_your_new_branch]

- Después de hacer cambios etc. a la branch, hacer un commit al repositorio local:

      $ git add

      $ git commit

- ["Subir" los cambios a GitHub](http://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/):

      $ git push origin [name_of_your_new_branch]

- Hacer un pull request con los cambios (desde la web).

- Volver a actualizar el repositorio local (master branch) y eliminar la rama (local):

$ git pull upstream master
$ git branch -d [name_of_your_new_branch]

- Actualizar de nuevo a `master` branch del repositorio local:

$ git push origin master
$ git push --delete origin [name_of_your_new_branch]

**Notas**:

$ git branch # para ver las branches
$ git checkout [branch] # cambiar a otra rama
