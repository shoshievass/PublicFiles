mutable struct Bid{T1,T2,T3,T4}
    init_index::Int64
    a::T1
    b::T2
    q::T3
    a_over_q::T4
end

function quadprog_fastSolve(a,b,q, s)
    num_bids = length(a)
    bidvec = [Bid(i, a[i],b[i],q[i],(a[i]/q[i])) for i=1:num_bids]

    sorted_bids = sort(bidvec, by = x -> x.a_over_q, rev=true)
    
    function x_star(bid, v)
        out = (bid.a + v * bid.q)/(2*bid.b) 
    end

    function v_k_star(k_plus_one)
        include = [i < k_plus_one ? 1 : 0 for i=1:num_bids]
        sum_1 = sum( include[i] * (sorted_bids[i].a * sorted_bids[i].q)/(sorted_bids[i].b)  for i=1:num_bids)
        sum_2 = sum( include[i] * (sorted_bids[i].q * sorted_bids[i].q)/(sorted_bids[i].b)  for i=1:num_bids)
        out_temp = sum(include) > 0 ? (2*s - sum_1)/sum_2 : 0
        out = k_plus_one <= num_bids ? min(out_temp, -sorted_bids[k_plus_one].a_over_q) : out_temp
        return(out)
    end

    function dual_eval(v)
        include = [-bid.a_over_q <= v ? 1 : 0 for bid in sorted_bids]
        x_stars = [include[i]*x_star(sorted_bids[i], v) for i=1:num_bids]
        out = sum( include[i] * ((sorted_bids[i].a + v*sorted_bids[i].q)^2/(4*sorted_bids[i].b) ) for i=1:num_bids) - v*s
        return(out, x_stars)
    end
        
    v_ks = [v_k_star(k) for k=1:(num_bids+1)]
    duals = [dual_eval(v_k) for v_k in v_ks]
    
    sort!(duals, by = x -> x[1])
    
    opt = duals[1]
    x_star_sorted = opt[2]
    x_star_orig_order = copy(x_star_sorted)
    for i = 1:num_bids
        bid = sorted_bids[i]
        x_star_orig_order[bid.init_index] = x_star_sorted[i]
    end
    
   return(x_star_orig_order)   
    
end

