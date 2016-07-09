using QuantEcon
using PyPlot

function get_simultaneous_P(epsilon, N, p)
    P = Array(Float64, (N, N))
    for i in 1:N
        if i < p * N
            P[i, 1:end] = [binomial(N, j)*(epsilon)^j*(1-epsilon)^(N-j) for j in 1:N]
            P[i, end] = 1-sum(P[i, 1:end-1])
        elseif i == p * N
            P[i, 1:end] = [binomial(N, j)*(1/2)^N for j in 1:N]
            P[i, end] = 1-sum(P[i, 1:end-1])
        else
            P[i, 1:end] = [binomial(N, j)*(epsilon)^(N-j)*(epsilon)^j for j in 1:N]
            P[i, end] = 1-sum(P[i, 1:end-1])
        end
    end
    return P
end

function get_sequential_P(epsilon, N, p)
    P = zeros(Float64, (N, N))
    for i in 1:N
        if i == 1
            P[i, i + 1] = (N - i)/N*epsilon/2
            if i > p * (N-1)
                P[i, i + 1] += 1 - epsilon
            elseif i == p * (N-1)
                P[i, i + 1] += 1/2
            end
            P[i, i] = 1 - P[i, i + 1]
        elseif i == N
            P[i, i-1] = i/N*epsilon/2
            if i - 1 == p * (N - 1)
                P[i, i-1] = 1/2
            elseif i - 1 < p * (N - 1)
                P[i, i-1] = 1 - epsilon
            end
            P[i, i] = 1 - P[i, i-1]
        else
            P[i, i+1] = (N - i)/N*epsilon/2
            if 1 > p * (N - i)
                P[i, i+1] += 1 - epsilon
            elseif 1 == p * (N - i)
                P[i, i+1] += 1/2
            end
            P[i, i-1] = i/N*epsilon/2
            if i - 1 == p * (N - 1)
                P[i, i-1] = 1/2
            elseif i - 1 < p * (N - 1)
                P[i, i-1] = 1 - epsilon
            end
            P[i, i] = 1 - P[i, i-1] - P[i, i + 1]
        end
    end
    return P
end

type KMR
    mc
    N
    epsilon
    Xs
end

get_p_2x2(M) = (M[1, 1, 1] - M[2, 1, 1])/(M[1, 1, 1] - M[1, 2, 1] - M[2, 1, 1] + M[2, 2, 1])
KMR(M, N::Int, epsilon::Float64; simul=true) = KMR(simul ? MarkovChain(get_simultaneous_P(epsilon, N, get_p_2x2(M))) : MarkovChain(get_sequential_P(epsilon, N, get_p_2x2(M))), N, epsilon, None)

function simulate!(kmr::KMR, T; init=0)
    kmr.Xs = simulate(kmr.mc, T)
end

plot_sample_path(kmr::KMR) = plot(kmr.Xs)#, label="epsilon = $kmr.epsion")

M = Array(Int, (2, 2, 2))
M[1, 1, :] = [4, 4]
M[1, 2, :] = [0, 3]
M[2, 1, :] = [3, 0]
M[2, 2, :] = [2, 2]

N = 20
T = 1000
epsilon = 1/5
p = 1/3

"""
for e in linspace(1/2, 0, 7)
    P = get_sequential_P(e, N, p)
    mc = MarkovChain(P)
    X = simulate(mc, T)
    plot(X, label="epsilon = $e")
end
legend()
"""

for e in linspace(1/2, 0, 7)
    kmr = KMR(M, N, e)
    simulate!(kmr, T)
    plot_sample_path(kmr)
end
#println(Xs)
#plot(Xs)
#println(P)
#println(mean(X .== 2))
