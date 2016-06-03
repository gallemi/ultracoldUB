

Bruno- Ubuntu laptop

-------------------
primera compilacion
-------------------

Para poder correr la primera version de spinrgs4.mb he hecho:

python2.7 springs4.py 

al principio me ha dado errores de blas, he reinstalado algunos paquetes de blas
sudo apt-get install libblas3
y he tenido que quitar algunos
sudo apt-get remove libopenblas-base

Con algun tira y afloja al final ha corrido perfectamente.


---------------------
putadillas al codigo
---------------------

1) Cambiando la masa a 5 kilos se descuadra todo
