f(x) = a*x+b
g(x) = c*x**2+d*x+e
fit f(x) "timesNewman.txt" using 1:2 via a, b
fit g(x) "timesNewman.txt" using 1:2 via c,d,e
plot f(x), g(x), "timesNewman.txt" using 1:2
