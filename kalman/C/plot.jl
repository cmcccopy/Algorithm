using CairoMakie

using DelimitedFiles

accy=readdlm("accy.txt");
mat=readdlm("matt.txt");

f = Figure()
aax = Axis(f[1, 1])
x=1:1:1713
lines!(aax,x,[accy...],color = :red)
lines!(aax,x,[mat...],color = :blue)

position=readdlm("position.txt")
apx = Axis(f[2, 1])
lines!(apx,x,[position...],color = :red)

speed=readdlm("speed.txt")
asx = Axis(f[3, 1])
lines!(asx,x,[speed...],color = :red)

save("1.png",f)