using LinearAlgebra
using Printf

D = [0.0313 2.5 7.0 5.6569;
    2.5 0.0313 4.5 4.272;
    7.0 4.5 0.0313 5.0;
    5.6569 4.272 5.0 0.0081]

Zabc = zeros(Complex, 4, 4)
for i in 1:4
    for j in 1:4
        if i != j
            Zabc[i, j] = 0.0953 + 0.12134 * (log(1.0 / D[i, j]) + 7.93402)im
        elseif i == 4
            Zabc[i, j] = 0.592 + 0.0953 + 0.12134 * (log(1.0 / D[i, j]) + 7.93402)im
        else
            Zabc[i, j] = 0.1859 + 0.0953 + 0.12134 * (log(1.0 / D[i, j]) + 7.93402)im
        end
    end
end

Zabc

Zabck = Zabc[1:3, 1:3] - Zabc[1:3, 4:4] * Zabc[4:4, 1:3] / Zabc[4, 4]

l = 10e3 / 5280.0

Zabck_ = Zabck * l

a = I
b = Zabck_
c = zeros(Complex, 3, 3)
d = a
A = inv(a)
B = A * b

p(m, a) = m * cis(deg2rad(a))

function dv(v)
    for i in v
        m = abs(i)
        a = rad2deg(angle(i))
        @printf "%.2f ∠%.2fº\n" m a
    end
end

va = 12.47e3 / √3

Vabc1_ = [p(va, 0.0); p(va, -120.0); p(va, 120.0);;]

dv(Vabc1_)

s = 2.0e6

Sabc = fill(p(s, acosd(0.9)), (3, 1))

begin
    Vabc2 = Vabc1_
    n = 0
    while n < 10
        global Vabc1, Iabc2, Vabc2, n
        n = n + 1
        Iabc2 = conj.(Sabc ./ Vabc2)

        # Backward Sweep
        Vabc1 = a * Vabc2 + b * Iabc2
        Iabc1 = c * Vabc2 + d * Iabc2

        error = maximum(abs.((Vabc1 - Vabc1_) ./ Vabc1_))
        @printf "Iteration %i absolute error = %.4f\n" n error
        if error < 1e-5
            break
        end

        # Forward Sweep
        Vabc1 = Vabc1_
        Vabc2 = A * Vabc1 - B * Iabc2
    end
end

dv(Vabc2)

dv(Iabc2)

dv(Vabc1)

abs(sum(Iabc2))

begin
    tn = -Zabc[4:4, 1:3] / Zabc[4, 4]
    In2 = tn * Iabc2
end

Ig2 = -[sum(Iabc2)] - In2

dv(In2 + Ig2)

