f(x) = a*x+b
g(x) = c*x**2+d*x+e
fit f(x) "timesRenamingAfter.txt" using 2:4 via a, b
fit g(x) "timesRenamingAfter.txt" using 2:4 via c,d,e
plot f(x), g(x), "timesRenamingAfter.txt" using 2:4
