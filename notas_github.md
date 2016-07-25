**Pasos para actualizar el repositorio local y hacer un pull request con los cambios:**
- [Sincronizar el repositorio local (fork) con el original](https://help.github.com/articles/syncing-a-fork/):

      $ git fetch upstream
      $ git checkout master
      $ git merge upstream/master

- [Actualizar mi repositorio con el local](https://help.github.com/articles/pushing-to-a-remote/):

      $ git push origin master

- Para actualizar el repositorio con los cambios locales:

      $ git status
      $ git add
      $ git commit
      $ git push origin master
      $ git status  # para verificar

- El pull request se puede hacer desde la web.

**Notas:**
- Para guardar y salir del editor vim: `Esc shift+Z shift+Z`.
- No tengo muy clara la utilidad de las branchs cuando tu repositorio es un fork de otro.
