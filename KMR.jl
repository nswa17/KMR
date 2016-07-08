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
        xs[i] = 3 - r
    else#
        xs[i] = r
    end
end

function simultaneous_revisions!(xs, N, M, epsilon)
    new_xs = copy(xs)
    for i in 1:length(xs)
        r = best_response(xs, i, N, M)
        if epsilon > rand()#突然変異
            new_xs[i] = 3 - r
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

M = Array(Int, (2, 2, 2))
M[1, 1, :] = [4, 4]
M[1, 2, :] = [0, 3]
M[2, 1, :] = [3, 0]
M[2, 2, :] = [2, 2]

N = 50
T = 100
epsilon = 0.05

xs = ones(Int, N)#######自由

Xs = main_sim(xs, M, N, T, epsilon)
println(Xs)
