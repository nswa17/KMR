using QuantEcon
using PyPlot

function get_simultaneous_P(epsilon, N, p)
    P = Array(Float64, (N+1, N+1))
    for i in 1:(N+1)
        if i - 1 < p * N
            P[i, 1:end] = [binomial(N, j)*(epsilon)^j*(1-epsilon)^(N-j) for j in 0:N]
            P[i, end] = 1-sum(P[i, 1:end-1])
        elseif i - 1 == p * N
            P[i, 1:end] = [binomial(N, j)*(1/2)^N for j in 0:N]
            P[i, end] = 1-sum(P[i, 1:end-1])
        else
            P[i, 1:end] = [binomial(N, j)*(epsilon)^(N-j)*(epsilon)^j for j in 0:N]
            P[i, end] = 1-sum(P[i, 1:end-1])
        end
    end
    return P
end

function get_sequential_P(epsilon, N, p)
    P = zeros(Float64, (N+1, N+1))
    for i in 1:N+1
        if i == 1
            P[i, i + 1] = (N - (i-1))/N*epsilon/2
            if (i-1) > p * (N-1)
                P[i, i + 1] += 1 - epsilon
            elseif (i-1) == p * (N-1)
                P[i, i + 1] += 1/2
            end
            P[i, i] = 1 - P[i, i + 1]
        elseif i == N+1
            P[i, i-1] = (i-1)/N*epsilon/2
            if i - 2 == p * (N - 1)
                P[i, i-1] = 1/2
            elseif i - 2 < p * (N - 1)
                P[i, i-1] = 1 - epsilon
            end
            P[i, i] = 1 - P[i, i-1]
        else
            P[i, i+1] = (N - (i-1))/N*epsilon/2
            if 1 > p * (N - (i-1))
                P[i, i+1] += 1 - epsilon
            elseif 1 == p * (N - (i-1))
                P[i, i+1] += 1/2
            end
            P[i, i-1] = (i-1)/N*epsilon/2
            if i - 2 == p * (N - 1)
                P[i, i-1] = 1/2
            elseif i - 2 < p * (N - 1)
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
    kmr.Xs = simulate(kmr.mc, T, init+1) .- 1
end

plot_sample_path(kmr::KMR) = plot(kmr.Xs)#, label="epsilon = $kmr.epsion")

function plot_stationary_dist(kmr::KMR)
    fig, ax = subplots()
    #ax[:hist](mc_compute_stationary(kmr.mc))
    ax[:bar](collect(1:kmr.N), mc_compute_stationary(kmr.mc))
end

function plot_empirical_dist(kmr::KMR)
    fig, ax = subplots()
    ax[:hist](kmr.Xs, normed=true)
end

function plot_sample_path(kmr::KMR)
    plot(kmr.Xs)
end

M = Array(Int, (2, 2, 2))
M[1, 1, :] = [4, 4]
M[1, 2, :] = [0, 3]
M[2, 1, :] = [3, 0]
M[2, 2, :] = [2, 2]

N = 20
T = 100
epsilon = 1/5

kmr = KMR(M, N, epsilon)
simulate!(kmr, T)
#plot_stationary_dist(kmr)
plot_empirical_dist(kmr)
#plot_sample_path(kmr)

"""
for e in linspace(1/2, 0, 7)
    kmr = KMR(M, N, e)
    simulate!(kmr, T)
    #plot_sample_path(kmr)
    #plot_empirical_dist(kmr)
    plot_stationary_dist(kmr)
end
#println(Xs)
#plot(Xs)
#println(P)
#println(mean(X .== 2))
"""
