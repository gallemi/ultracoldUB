reset
set term gif size 1200,600 animate  delay 30 loop 1 optimize 
set output "minus.gif"
i=1
load "gif_minus.gnu"
