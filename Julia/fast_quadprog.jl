using ForwardDiff

mutable struct Bid{T1,T2,T3,T4}
    init_index::Int64
    a::T1
    b::T2
    q::T3
    a_over_q::T4
end

function quadprog_exactSolve(a,b,q, s, detail = false)
    num_bids = length(a)
    bidvec = [Bid(i, a[i],b[i],q[i],(a[i]/q[i])) for i=1:num_bids ]
    sorted_bids = sort(bidvec, by = x -> x.a_over_q, rev=true)
    
    if detail
        println("sorted_bids: ")
        for bid in sorted_bids
            println(bid)
        end
        
    end
    
    ## Save sorted vector for faster calculations
    sorted_a = [sorted_bids[i].a for i=1:num_bids]
    sorted_b = [sorted_bids[i].b for i=1:num_bids]
    sorted_q = [sorted_bids[i].q for i=1:num_bids]
    min_sorted_a_over_q = [-sorted_bids[i].a_over_q for i=1:num_bids]
    min_sorted_a_over_q_nxt = map(x->zero(float(x)), sorted_a)
#     @show min_sorted_a_over_q_nxt = [0.0 for i = 1:num_bids]
    
    for i=1:(num_bids-1)
        min_sorted_a_over_q_nxt[i] = min_sorted_a_over_q[i+1]
    end
    min_sorted_a_over_q_nxt[num_bids] = Inf


    ## Compute v_ks
    t1 = sorted_a .* sorted_q ./ sorted_b
    t2 = sorted_q.^2 ./ sorted_b
    t1_acc = accumulate(+, t1)
    t2_acc = accumulate(+, t2)
    
    ## Set v_ks and handle cases where b_i = 0
    xst_try = (2*s .- t1_acc) ./ t2_acc
    v_ks = copy(xst_try)
    v_ks = replace!(v_ks, NaN => Inf) 
    
    
    for i=1:(num_bids)
        v_ks[i] = sorted_b[i] <= 0 ? min_sorted_a_over_q[i] : min.(xst_try[i], min_sorted_a_over_q_nxt[i])
        if sorted_b[i] <= 0 ## only want to set the _first_ non-quad b
            break
        end
    end
    
    ## Compute duals
    function dual_eval(v)
#         x_stars = [0.0 for  i=1:num_bids]
#         x_stars = map(x->zero(float(x)), sorted_a)
        x_stars = map(x->zero(float(eltype(v))), 1:num_bids)
        dual_eval = 0.
        cum_qx = 0.
        for i=1:num_bids
            if min_sorted_a_over_q[i] > v
                break
            end
            if sorted_b[i] <= 0.
                x_stars[i] = (s-cum_qx)/sorted_q[i] 
                dual_eval += ( (sorted_a[i]/sorted_q[i]) + v ) * (s - cum_qx)
                break
            else
                x_stars[i] = ( (sorted_a[i] + (v * sorted_q[i]) ) / (2 * sorted_b[i]) )
                dual_eval += ( (sorted_a[i] + (v * sorted_q[i]) )^2 / (4 * sorted_b[i]) )
                cum_qx +=  sorted_q[i]*x_stars[i]
            end
        end

        out = dual_eval - v*s
        return(out, x_stars)
    end
    
    duals = [dual_eval(v_k) for v_k in v_ks]
    
    
    if detail
        println("v_ks: ", v_ks)
        println("dual obj: ", [duals[i][1] for i=1:length(duals)])
    end
    
    ## check corner solutions
    function getCornerObj(i)
        x_i_st = s / sorted_q[i]
        obj = sorted_a[i]*x_i_st - sorted_b[i]*(x_i_st^2)
        return((obj,i))
    end
    
    corner_vals = [getCornerObj(i) for i=1:num_bids]

    
    sort!(duals, by = x -> x[1])
    sort!(corner_vals, by = x -> x[1], rev=true)

    duals_opt = duals[1]
    corners_opt = corner_vals[1]
    
    if duals_opt[1] >= corners_opt[1]
        opt_val = duals_opt[1]
        x_star_sorted = duals_opt[2]
    else
        opt_val = corners_opt[1]
        opt_i = corners_opt[2]
        x_star_sorted = map(x->zero(float(x)), sorted_a)
        x_star_sorted[opt_i] = (s / sorted_q[opt_i])
    end
    x_star_orig_order = copy(x_star_sorted)
    for i = 1:num_bids
        bid = sorted_bids[i]
        x_star_orig_order[bid.init_index] = x_star_sorted[i]
    end
    
    return(x_star_orig_order)
    
end