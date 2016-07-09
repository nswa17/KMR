using PyPlot

function get_X(xs)
    return sum(filter(x->x==2, xs))
end

function random_choice(xs)
    return rand(collect(1:length(xs)))
end

function best_response(xs, i, N, M)
    if xs[i] == 1
        p = get_X(xs)/(N-1)
    else
        p = (get_X(xs) - 1)/(N-1)
    end
    if p * M[1, 2, 1] + (1-p) * M[1, 1, 1] > p * M[2, 2, 1] + (1-p) * M[2, 1, 1]
        return 1
    elseif p * M[1, 2, 1] + (1-p) * M[1, 1, 1] < p * M[2, 2, 1] + (1-p) * M[2, 1, 1]
        return 2
    else
        return rand([1, 2])
    end
end

function sequential_revisions!(xs, N, M, epsilon)
    i = rand(collect(1:length(xs)))
    r = best_response(xs, i, N, M)
    if epsilon > rand()#突然変異
        if 1/2 > rand()
            xs[i] = 3 - r
        else
            xs[i] = r
        end
    else#
        xs[i] = r
    end
end

function simultaneous_revisions!(xs, N, M, epsilon)
    new_xs = copy(xs)
    for i in 1:length(xs)
        r = best_response(xs, i, N, M)
        if epsilon > rand()#突然変異
            if 1/2 > rand()
                new_xs[i] = 3 - r
            else
                new_xs[i] = r
            end
        else#
            new_xs[i] = r
        end
    end
    xs = new_xs
end

function main_sim(xs, M, N, T, epsilon; simul = false)
    Xs = Array(Int, T)
    for t in 1:T
        if simul
            simultaneous_revisions!(xs, N, M, epsilon)
        else
            sequential_revisions!(xs, N, M, epsilon)
        end
        Xs[t] = get_X(xs)
        #println(xs)
    end
    return Xs
end

type KMR
    xs
    M
    N
    T
    epsilon
    Xs
end

KMR(M, N::Int, epsilon) = KMR(ones(Int, N), M, N, 0, epsilon, Array(Int, T))

function simulate!(kmr::KMR, T; init = 0)
    kmr.xs[1:init] = 2
    kmr.Xs = main_sim(kmr.xs, M, N, T, epsilon)
end

function plot_sample_path(kmr)
    plot(kmr.Xs)
end

function plot_empirical_dist(kmr)
    fig, ax = subplots()
    ax[:hist](kmr.Xs)
end

M = Array(Int, (2, 2, 2))
M[1, 1, :] = [4, 4]
M[1, 2, :] = [0, 3]
M[2, 1, :] = [3, 0]
M[2, 2, :] = [2, 2]

N = 20
T = 1000
epsilon = 0.06
"""
loop = 10
for i in 1:loop
    xs = ones(Int, N)#######自由
    Xs = main_sim(xs, M, N, T, epsilon)
    plot(Xs)
end
"""
kmr = KMR(M, N, epsilon)
simulate!(kmr, T)
plot_sample_path(kmr)
plot_empirical_dist(kmr)
