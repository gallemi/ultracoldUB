reset
set term gif size 1000,800 animate
set term gif animate #optimize
set output "bright_soliton.gif"
i=1
load "animation.gnu"